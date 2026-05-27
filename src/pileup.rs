//! Streaming position-ordered pileup over raw BAM bytes.
//!
//! We decode each BAM record's payload directly (same approach as bam-head and
//! bam-depad) to avoid noodles' high-level record API, which has a different
//! surface area per version. Each read is expanded into an `ActiveRead` that
//! maps reference positions to query base+quality observations, then we advance
//! position by position producing per-site consensus values.
//!
//! ## BAM payload layout (after 4-byte block_size)
//! ```text
//! Offset  Field
//!  0..4   ref_id     (i32, LE)
//!  4..8   pos        (i32, LE, 0-based)
//!  8      l_read_name (u8)
//!  9      mapq       (u8)
//! 10..12  bin        (u16 LE)
//! 12..14  n_cigar    (u16 LE)
//! 14..16  flag       (u16 LE)
//! 16..20  l_seq      (u32 LE)
//! 20..24  next_ref   (i32 LE)
//! 24..28  next_pos   (i32 LE)
//! 28..32  tlen       (i32 LE)
//! 32..    name (l_read_name bytes incl NUL)
//!         cigar (n_cigar * 4 bytes)
//!         seq   (ceil(l_seq/2) bytes, 4-bit packed)
//!         qual  (l_seq bytes)
//!         aux   ...
//! ```

use std::io::Read;

use rsomics_bamio::raw::{self, RawRecord};
use rsomics_common::Result;

use crate::errmod::{self, Errmod};

// BAM payload field offsets.
const REF_ID: usize = 0;
const POS: usize = 4;
const L_READ_NAME: usize = 8;
const MAPQ: usize = 9;
const N_CIGAR: usize = 12;
const FLAG: usize = 14;
const L_SEQ: usize = 16;
const FIXED_HEAD: usize = 32;

// BAM flag bits to skip (unmapped | secondary | qcfail | dup).
const SKIP_FLAGS: u16 = 0x4 | 0x100 | 0x200 | 0x400;
// Reverse-complement flag.
const FLAG_REV: u16 = 0x10;

// CIGAR op codes.
const OP_M: u8 = 0; // match/mismatch
const OP_I: u8 = 1; // insertion to ref
const OP_D: u8 = 2; // deletion from ref
const OP_N: u8 = 3; // skip
const OP_S: u8 = 4; // soft clip
const OP_EQ: u8 = 7; // seq match
const OP_X: u8 = 8; // seq mismatch

/// BAM 4-bit seq nibble → ACGT index (0-3, 255=ambiguous).
/// Table maps nibble values 0..15 (=ACMGRSVTWYHKDBN) → 0=A,1=C,2=G,3=T,255=ambig.
#[rustfmt::skip]
const NT4: [u8; 16] = [
    255, // 0 = '='
    0,   // 1 = A
    1,   // 2 = C
    255, // 3 = M
    2,   // 4 = G
    255, // 5 = R
    255, // 6 = S
    255, // 7 = V
    3,   // 8 = T
    255, // 9 = W
    255, // 10 = Y
    255, // 11 = H
    255, // 12 = K
    255, // 13 = D
    255, // 14 = B
    255, // 15 = N
];

/// A decoded active read in the pileup window.
struct ActiveRead {
    tid: i32,
    ref_start: i32, // 0-based
    is_rev: bool,
    mapq: u8,
    /// For each reference position offset from ref_start:
    ///   `query_idx[i]` = index into `bases`/`quals`, or `i32::MIN` for del/refskip.
    query_idx: Vec<i32>,
    /// Decoded bases: 0=A,1=C,2=G,3=T,255=ambiguous.
    bases: Vec<u8>,
    /// Per-base quality scores (raw phred, unclamped).
    quals: Vec<u8>,
}

impl ActiveRead {
    fn ref_end(&self) -> i32 {
        self.ref_start + self.query_idx.len() as i32
    }

    /// Retrieve the errmod base encoding for reference position `ref_pos`,
    /// or `None` if deleted/refskip/ambiguous/below min_baseQ.
    fn base_enc(&self, ref_pos: i32, min_base_q: u8) -> Option<u16> {
        let offset = (ref_pos - self.ref_start) as usize;
        let qi = *self.query_idx.get(offset)?;
        if qi < 0 {
            return None;
        }
        let qi = qi as usize;
        let b = *self.bases.get(qi)?;
        if b > 3 {
            return None;
        }
        let bq = *self.quals.get(qi)?;
        if bq < min_base_q {
            return None;
        }
        // Quality clamped to [4, 63] for errmod; also capped by mapq.
        let q = bq.min(self.mapq).clamp(4, 63);
        let strand = u16::from(self.is_rev);
        Some((q as u16) << 5 | strand << 4 | b as u16)
    }

    /// Parse from raw BAM payload bytes (no block_size prefix).
    fn from_bytes(rec: &[u8]) -> Option<Self> {
        if rec.len() < FIXED_HEAD {
            return None;
        }
        let flags = u16::from_le_bytes([rec[FLAG], rec[FLAG + 1]]);
        if flags & SKIP_FLAGS != 0 {
            return None;
        }
        let tid = i32::from_le_bytes(rec[REF_ID..REF_ID + 4].try_into().ok()?);
        if tid < 0 {
            return None; // unmapped
        }
        let pos = i32::from_le_bytes(rec[POS..POS + 4].try_into().ok()?);
        let mapq = rec[MAPQ].clamp(4, 63);
        let n_cigar = u16::from_le_bytes([rec[N_CIGAR], rec[N_CIGAR + 1]]) as usize;
        let l_seq = u32::from_le_bytes(rec[L_SEQ..L_SEQ + 4].try_into().ok()?) as usize;
        let l_name = rec[L_READ_NAME] as usize;
        let is_rev = flags & FLAG_REV != 0;

        let cigar_start = FIXED_HEAD + l_name;
        let seq_start = cigar_start + n_cigar * 4;
        let qual_start = seq_start + l_seq.div_ceil(2);

        if rec.len() < qual_start + l_seq {
            return None;
        }

        // Decode seq.
        let mut bases = Vec::with_capacity(l_seq);
        for i in 0..l_seq {
            let byte = rec[seq_start + i / 2];
            let nibble = if i % 2 == 0 { byte >> 4 } else { byte & 0x0f };
            bases.push(NT4[nibble as usize]);
        }

        // Decode qual.
        let quals = rec[qual_start..qual_start + l_seq].to_vec();

        // Build ref→query map from CIGAR.
        let mut query_idx: Vec<i32> = Vec::new();
        let mut q: i32 = 0;
        for ci in 0..n_cigar {
            let raw = u32::from_le_bytes(
                rec[cigar_start + ci * 4..cigar_start + ci * 4 + 4]
                    .try_into()
                    .ok()?,
            );
            let op = (raw & 0xf) as u8;
            let ol = (raw >> 4) as i32;
            match op {
                OP_M | OP_EQ | OP_X => {
                    for _ in 0..ol {
                        query_idx.push(q);
                        q += 1;
                    }
                }
                OP_D | OP_N => {
                    for _ in 0..ol {
                        query_idx.push(i32::MIN); // deletion / refskip
                    }
                }
                OP_I | OP_S => {
                    q += ol; // consumes query, not reference
                }
                _ => {} // H, P: skip
            }
        }

        Some(Self {
            tid,
            ref_start: pos,
            is_rev,
            mapq,
            query_idx,
            bases,
            quals,
        })
    }
}

/// Result of one `PileupEngine::step` call.
pub enum StepResult {
    /// A pileup position was computed.
    PileupPosition { tid: i32, pos: i32, site_cns: u16 },
    /// End of input.
    Done,
}

/// Streaming pileup engine — call `step()` repeatedly until `Done`.
pub struct PileupEngine<'a> {
    min_base_q: u8,
    em: &'a Errmod,
    active: Vec<ActiveRead>,
    cur_tid: i32,
    cur_pos: i32,
    /// Reads read ahead but not yet activated (sorted by tid, ref_start).
    pending: Vec<ActiveRead>,
    exhausted: bool,
}

impl<'a> PileupEngine<'a> {
    pub fn new(min_base_q: u8, em: &'a Errmod) -> Self {
        Self {
            min_base_q,
            em,
            active: Vec::new(),
            cur_tid: -1,
            cur_pos: 0,
            pending: Vec::new(),
            exhausted: false,
        }
    }

    /// Feed one raw record payload (no block_size prefix). Returns true if the
    /// record was consumed into the active set; false if it was queued as pending.
    fn ingest(&mut self, bytes: &[u8]) -> bool {
        let Some(ar) = ActiveRead::from_bytes(bytes) else {
            return true; // skip (unmapped/filtered)
        };
        if ar.tid < self.cur_tid || (ar.tid == self.cur_tid && ar.ref_start <= self.cur_pos) {
            self.active.push(ar);
            true
        } else {
            // Read is ahead of current position — park in pending.
            self.pending.push(ar);
            false
        }
    }

    /// Load all reads from `reader` whose start ≤ `cur_pos` into `active`.
    fn fill_active<R: Read>(&mut self, rec: &mut RawRecord, reader: &mut R) -> Result<()> {
        // First drain pending reads that now apply.
        let mut i = 0;
        while i < self.pending.len() {
            let p = &self.pending[i];
            if p.tid < self.cur_tid || (p.tid == self.cur_tid && p.ref_start <= self.cur_pos) {
                let ar = self.pending.remove(i);
                self.active.push(ar);
            } else {
                i += 1;
            }
        }

        // Then read new records until we hit one that starts after cur_pos.
        if self.exhausted {
            return Ok(());
        }
        loop {
            let n = raw::read_record(reader, rec)?;
            if n == 0 {
                self.exhausted = true;
                break;
            }
            let consumed = self.ingest(rec.as_bytes());
            if !consumed {
                break;
            }
        }
        Ok(())
    }

    /// Advance the engine to the first covered position.  Must be called once
    /// before the first `step` when we don't yet know `cur_tid`/`cur_pos`.
    fn bootstrap<R: Read>(&mut self, rec: &mut RawRecord, reader: &mut R) -> Result<bool> {
        if self.exhausted {
            return Ok(false);
        }
        loop {
            let n = raw::read_record(reader, rec)?;
            if n == 0 {
                self.exhausted = true;
                return Ok(false);
            }
            let bytes = rec.as_bytes();
            if bytes.len() < FIXED_HEAD {
                continue;
            }
            let flags = u16::from_le_bytes([bytes[FLAG], bytes[FLAG + 1]]);
            if flags & SKIP_FLAGS != 0 {
                continue;
            }
            let tid = i32::from_le_bytes(bytes[REF_ID..REF_ID + 4].try_into().unwrap());
            if tid < 0 {
                continue;
            }
            let pos = i32::from_le_bytes(bytes[POS..POS + 4].try_into().unwrap());
            self.cur_tid = tid;
            self.cur_pos = pos;
            self.ingest(bytes);
            return Ok(true);
        }
    }

    /// Step to the next covered position and return its consensus.
    /// On the very first call, reads until the first covered base.
    pub fn step<R: Read>(&mut self, rec: &mut RawRecord, reader: &mut R) -> Result<StepResult> {
        // Bootstrap on first call.
        if self.cur_tid < 0 && !self.bootstrap(rec, reader)? {
            return Ok(StepResult::Done);
        }

        loop {
            // Purge reads that end at or before cur_pos.
            self.active.retain(|r| r.ref_end() > self.cur_pos);

            // Load reads that start ≤ cur_pos.
            self.fill_active(rec, reader)?;

            // Check if any active read covers cur_pos.
            let covered = self.active.iter().any(|r| {
                r.tid == self.cur_tid && r.ref_start <= self.cur_pos && r.ref_end() > self.cur_pos
            });

            if !covered {
                // No coverage at cur_pos. Jump to the next covered position.
                let next = self.next_covered_position();
                match next {
                    Some((ntid, npos)) => {
                        self.cur_tid = ntid;
                        self.cur_pos = npos;
                        continue;
                    }
                    None => return Ok(StepResult::Done),
                }
            }

            // Build base observations.
            let mut bases: Vec<u16> = Vec::new();
            for ar in &self.active {
                if ar.tid == self.cur_tid
                    && ar.ref_start <= self.cur_pos
                    && ar.ref_end() > self.cur_pos
                    && let Some(enc) = ar.base_enc(self.cur_pos, self.min_base_q)
                {
                    bases.push(enc);
                }
            }

            let site_cns = if bases.is_empty() {
                0u16
            } else {
                errmod::gencns(self.em, &mut bases)
            };

            let result = StepResult::PileupPosition {
                tid: self.cur_tid,
                pos: self.cur_pos,
                site_cns,
            };

            // Advance to next position: min of (cur_pos+1 if any read extends, or next pending start).
            let max_end = self
                .active
                .iter()
                .filter(|r| r.tid == self.cur_tid)
                .map(|r| r.ref_end())
                .max();

            match max_end {
                Some(end) if end > self.cur_pos + 1 => {
                    self.cur_pos += 1;
                }
                _ => {
                    // Advance past end of current reads.
                    match self.next_covered_position() {
                        Some((ntid, npos)) => {
                            self.cur_tid = ntid;
                            self.cur_pos = npos;
                        }
                        None => {
                            // Signal done after returning this last position.
                            self.cur_tid = -1;
                        }
                    }
                }
            }

            return Ok(result);
        }
    }

    /// Find the smallest (tid, pos) covered by any active or pending read
    /// that is strictly after `(cur_tid, cur_pos)`.
    fn next_covered_position(&self) -> Option<(i32, i32)> {
        let cur = (self.cur_tid, self.cur_pos);

        // Candidates from active reads (their first reference position ≥ current).
        let from_active = self
            .active
            .iter()
            .filter(|r| {
                let rstart = (r.tid, r.ref_start);
                rstart > cur || (r.tid == self.cur_tid && r.ref_end() > self.cur_pos + 1)
            })
            .map(|r| {
                if r.tid == self.cur_tid && r.ref_end() > self.cur_pos + 1 {
                    (r.tid, self.cur_pos + 1)
                } else {
                    (r.tid, r.ref_start)
                }
            });

        // Candidates from pending reads.
        let from_pending = self
            .pending
            .iter()
            .filter(|r| (r.tid, r.ref_start) > cur)
            .map(|r| (r.tid, r.ref_start));

        from_active.chain(from_pending).min()
    }
}

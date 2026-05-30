//! `samtools targetcut` port: identify amplicon target intervals from pileup depth
//! and emit them as SAM records.
//!
//! ## Algorithm (mirrors `cut_target.c`)
//!
//! 1. Stream the BAM through a pileup engine, skipping unmapped/secondary/qcfail/dup reads.
//! 2. At each reference position compute a per-position consensus value (`u16`) via
//!    the MAQ error model (`errmod`): encodes depth (high byte, 0=no coverage) and
//!    the best-base quality/identity (low byte).
//! 3. After processing a chromosome, run a 2-state Viterbi HMM over the consensus
//!    array to label each position 0 (off-target) or 1 (on-target).
//! 4. Each run of consecutive on-target positions is emitted as a single SAM record:
//!    QNAME=`chr:start-end` (1-based), FLAG=0, RNAME=chr, POS=start (1-based),
//!    MAPQ=60, CIGAR=`NM`, SEQ=consensus bases (ACGT or N), QUAL=phred chars.
//!
//! The reference FASTA (`-f`) is optional; when supplied it enables BAQ
//! (`sam_prob_realn`) which samtools calls before the pileup. We do not
//! implement BAQ in this port (it requires htslib internals); the `-f` flag is
//! accepted for compatibility but silently ignored with a warning.

mod errmod;
mod pileup;

use std::io::{BufWriter, Write};
use std::num::NonZero;
use std::path::Path;

use noodles::sam;
use rsomics_bamio::raw::RawRecord;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

// HMM parameter defaults from cut_target.c `g_param`.
// e[state][coverage_class]: emission score for being in state `state` observing class c.
// p[from][to]: transition score (p[0][1] = cost of entering the target, default -14000).
// coverage class c: 0 = no coverage (cns==0), 1 = low depth (cns>>8 == 0, i.e. depth 1..255
//   but the depth byte is 0, meaning depth<1 — actually depth field is the second byte;
//   class 2 = high depth (cns>>8 > 0, i.e. depth in [1..255]).
//
// From the source:
//   e[0] = {0,0,0};     // off-target emission for classes 0,1,2
//   e[1] = {-4,1,6};    // on-target emission for classes 0,1,2
//   p[0][0]=0, p[0][1]=-14000  (0→0, 0→1)
//   p[1][0]=0, p[1][1]=0       (1→0, 1→1)

#[derive(Debug, Clone)]
pub struct HmmParams {
    /// e[state][coverage_class]
    pub e: [[i32; 3]; 2],
    /// p[from][to]
    pub p: [[i32; 2]; 2],
}

impl Default for HmmParams {
    fn default() -> Self {
        Self {
            e: [[0, 0, 0], [-4, 1, 6]],
            p: [[0, -14000], [0, 0]],
        }
    }
}

#[derive(Debug, Clone)]
pub struct TargetcutOpts {
    pub min_base_q: u8,
    pub in_penalty: i32,
    pub em0: Option<i32>,
    pub em1: Option<i32>,
    pub em2: Option<i32>,
    /// Accepted for CLI compat; BAQ re-alignment is not implemented (see module doc).
    pub reference: Option<std::path::PathBuf>,
}

impl Default for TargetcutOpts {
    fn default() -> Self {
        Self {
            min_base_q: 13,
            in_penalty: 14000,
            em0: None,
            em1: None,
            em2: None,
            reference: None,
        }
    }
}

#[derive(Debug, Default, Clone, Serialize)]
pub struct TargetcutStats {
    pub targets: u64,
}

pub fn targetcut(
    input: &Path,
    output_path: Option<&Path>,
    opts: &TargetcutOpts,
    workers: NonZero<usize>,
) -> Result<TargetcutStats> {
    match output_path {
        Some(path) => {
            let file = std::fs::File::create(path).map_err(|e| {
                RsomicsError::InvalidInput(format!("creating {}: {e}", path.display()))
            })?;
            let mut out = BufWriter::with_capacity(256 * 1024, file);
            let stats = run(input, &mut out, opts, workers)?;
            out.flush().map_err(RsomicsError::Io)?;
            Ok(stats)
        }
        None => {
            let stdout = std::io::stdout();
            let mut out = BufWriter::with_capacity(256 * 1024, stdout.lock());
            let stats = run(input, &mut out, opts, workers)?;
            out.flush().map_err(RsomicsError::Io)?;
            Ok(stats)
        }
    }
}

fn run<W: Write>(
    input: &Path,
    out: &mut W,
    opts: &TargetcutOpts,
    workers: NonZero<usize>,
) -> Result<TargetcutStats> {
    if opts.reference.is_some() {
        eprintln!(
            "[targetcut] WARNING: -f/--reference accepted but BAQ re-alignment is not \
             implemented in this port; reference is ignored"
        );
    }

    let mut params = HmmParams::default();
    params.p[0][1] = -opts.in_penalty;
    if let Some(v) = opts.em0 {
        params.e[1][0] = v;
    }
    if let Some(v) = opts.em1 {
        params.e[1][1] = v;
    }
    if let Some(v) = opts.em2 {
        params.e[1][2] = v;
    }

    let em = errmod::Errmod::new(1.0 - errmod::ERR_DEP);

    let mut reader = rsomics_bamio::open_with_workers(input, workers)?;
    let header: sam::Header = reader.read_header().map_err(RsomicsError::Io)?;

    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
        .collect();

    let ref_lengths: Vec<usize> = header
        .reference_sequences()
        .values()
        .map(|sq| usize::from(sq.length()))
        .collect();

    let mut engine = pileup::PileupEngine::new(opts.min_base_q, &em);
    let mut stats = TargetcutStats::default();

    // Pileup is inherently single-threaded (position-ordered scan), so we use
    // the parallel reader's inner BGZF stream directly with a single RawRecord.
    // `get_mut()` returns &mut Box<dyn BufRead + Send> which implements Read.
    let inner = reader.get_mut();
    let mut raw_rec = RawRecord::default();

    let mut prev_tid: i32 = -1;
    // Per-chromosome consensus array; resized as chromosomes appear.
    let mut cns: Vec<u16> = Vec::new();

    while let pileup::StepResult::PileupPosition { tid, pos, site_cns } =
        engine.step(&mut raw_rec, inner)?
    {
        if tid != prev_tid {
            // Chromosome boundary: run HMM on completed chromosome.
            if prev_tid >= 0 && !cns.is_empty() {
                let chrom = ref_names
                    .get(prev_tid as usize)
                    .map(|s| s.as_str())
                    .unwrap_or("?");
                let n = process_cns(out, chrom, &cns, &params)?;
                stats.targets += n;
            }
            let new_len = ref_lengths.get(tid as usize).copied().unwrap_or(0);
            cns.clear();
            cns.resize(new_len, 0u16);
            prev_tid = tid;
        }
        if (pos as usize) < cns.len() {
            cns[pos as usize] = site_cns;
        }
    }

    // Flush last chromosome.
    if prev_tid >= 0 && !cns.is_empty() {
        let chrom = ref_names
            .get(prev_tid as usize)
            .map(|s| s.as_str())
            .unwrap_or("?");
        let n = process_cns(out, chrom, &cns, &params)?;
        stats.targets += n;
    }

    Ok(stats)
}

/// Run the 2-state Viterbi HMM over `cns[0..l]` and emit one SAM record
/// per detected on-target run. Returns the number of target records emitted.
///
/// Mirrors `process_cns` in `cut_target.c` exactly:
/// - Coverage class: 0 if cns==0, 1 if cns>>8==0 (depth byte 0 but base present), 2 otherwise.
/// - Forward fill, backtrack, then scan for on-target runs.
/// - Each on-target run [s, i) is emitted as a SAM record with:
///   QNAME=`chr:s+1-i`, FLAG=0, RNAME=chr, POS=s+1, MAPQ=60, CIGAR=`(i-s)M`,
///   SEQ=consensus ACGT/N chars, QUAL=phred chars.
fn process_cns<W: Write>(out: &mut W, chrom: &str, cns: &[u16], params: &HmmParams) -> Result<u64> {
    let l = cns.len();
    if l == 0 {
        return Ok(0);
    }

    // b[i]: bits[1:0]=best-from for state-0, bits[3:2]=best-from for state-1, bits[5:4]=final label.
    let mut b: Vec<u8> = vec![0u8; l];
    let mut f = [0i32; 2];

    for i in 0..l {
        let c = coverage_class(cns[i]);

        let t0 = f[0] + params.e[0][c] + params.p[0][0];
        let t1 = f[1] + params.e[0][c] + params.p[1][0];
        let (f0_new, b0) = if t0 >= t1 { (t0, 0u8) } else { (t1, 1u8) };

        let u0 = f[0] + params.e[1][c] + params.p[0][1];
        let u1 = f[1] + params.e[1][c] + params.p[1][1];
        let (f1_new, b1) = if u0 >= u1 { (u0, 0u8) } else { (u1, 1u8) };

        b[i] = b0 | (b1 << 1);
        f[0] = f0_new;
        f[1] = f1_new;
    }

    let mut s = if f[0] >= f[1] { 0u8 } else { 1u8 };
    for i in (1..l).rev() {
        b[i] |= s << 2; // store final label in bits[5:4]
        s = (b[i] >> s) & 1;
    }

    let mut start: Option<usize> = None;
    let mut targets = 0u64;

    #[allow(clippy::needless_range_loop)]
    for i in 0..=l {
        let label = if i < l { (b[i] >> 2) & 3 } else { 0 };

        if label == 0 {
            // Off-target position (or past end).
            if let Some(s) = start.take() {
                emit_target(out, chrom, s, i, cns)?;
                targets += 1;
            }
        } else if start.is_none() {
            start = Some(i);
        }
    }

    Ok(targets)
}

/// Coverage class for the Viterbi emission: 0=no coverage, 1=covered but depth byte 0, 2=depth byte >0.
/// Mirrors `main_cut_target` in `cut_target.c`.
#[inline]
fn coverage_class(cns: u16) -> usize {
    if cns == 0 {
        0
    } else if cns >> 8 == 0 {
        1
    } else {
        2
    }
}

/// Emit one SAM target-interval record to `out`.
/// Format matches `process_cns` in `cut_target.c`:
///   QNAME  chr:start+1-end   (1-based, inclusive)
///   FLAG   0
///   RNAME  chr
///   POS    start+1            (1-based)
///   MAPQ   60
///   CIGAR  (end-start)M
///   RNEXT  *
///   PNEXT  0
///   TLEN   0
///   SEQ    consensus ACGT chars (N for no coverage or depth=0)
///   QUAL   phred-scale quals from cns upper nibble
fn emit_target<W: Write>(
    out: &mut W,
    chrom: &str,
    start: usize,
    end: usize,
    cns: &[u16],
) -> Result<()> {
    let len = end - start;
    write!(
        out,
        "{}:{}-{}\t0\t{}\t{}\t60\t{}M\t*\t0\t0\t",
        chrom,
        start + 1,
        end,
        chrom,
        start + 1,
        len
    )
    .map_err(RsomicsError::Io)?;
    for &v in &cns[start..end] {
        let c = v >> 8;
        let base = if c == 0 {
            b'N'
        } else {
            b"ACGT"[(c & 3) as usize]
        };
        out.write_all(&[base]).map_err(RsomicsError::Io)?;
    }
    out.write_all(b"\t").map_err(RsomicsError::Io)?;
    for &v in &cns[start..end] {
        let q = (v >> 8 >> 2) as u8;
        out.write_all(&[33 + q]).map_err(RsomicsError::Io)?;
    }
    out.write_all(b"\n").map_err(RsomicsError::Io)?;
    Ok(())
}

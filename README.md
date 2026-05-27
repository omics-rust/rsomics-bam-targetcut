# rsomics-bam-targetcut

Identify amplicon target intervals from pileup depth — Rust port of
`samtools targetcut`.

```sh
rsomics-bam-targetcut amplicon.bam          # default HMM parameters
rsomics-bam-targetcut -Q 20 amplicon.bam    # raise minimum base quality
rsomics-bam-targetcut -i 100 amplicon.bam   # lower in-penalty for short amplicons
```

## What it does

1. Streams the BAM (skipping unmapped, secondary, QC-fail, and duplicate reads).
2. Builds a per-position consensus base + quality estimate via the MAQ error
   model (a pure Rust port of `htslib/errmod.c`).
3. Runs a 2-state Viterbi HMM over the consensus array to label each position
   as on-target (1) or off-target (0).
4. Emits one SAM record per detected on-target run:
   `QNAME=chr:start-end`, `FLAG=0`, `MAPQ=60`, `CIGAR=NM`, `SEQ=consensus`.

Output format is SAM text on stdout, byte-identical to `samtools targetcut`
on the same input.

## Options

| Flag | Meaning | Default |
|---|---|---|
| `-Q, --min-baseq INT` | Minimum base quality for pileup | 13 |
| `-i, --in-penalty INT` | HMM cost for entering a target interval | 14000 |
| `--em0 INT` | On-target emission score at zero coverage | -4 |
| `--em1 INT` | On-target emission score at low coverage | 1 |
| `--em2 INT` | On-target emission score at high coverage | 6 |
| `-f, --reference FILE` | Reference FASTA (accepted; BAQ not implemented) | — |
| `-o, --output FILE` | Output SAM (default stdout) | — |
| `-t, --threads INT` | BGZF decoder threads | 1 |

**Penalty note**: the default in-penalty (14000) requires a sufficiently long
amplicon at high depth to overcome the entry cost. For shorter amplicons
(< 500 bp), lower the penalty with `-i 100` to `-i 1000`.

## How it is fast

- Reads are decoded from raw BAM bytes (no noodles high-level decode path)
  into a compact `ref_start + query_idx + bases + quals` struct.
- The pileup engine is a streaming position-ordered scanner with an active
  read list; it never buffers the full chromosome.
- The MAQ error model is table-driven and computed once at startup;
  per-position consensus is a fixed-size array multiply-add, not a heap
  allocation.

## Origin

This crate is an independent Rust reimplementation of `samtools targetcut`
based on:
- The upstream MIT-licensed source (`cut_target.c`, `htslib/errmod.c`):
  the HMM parameter defaults, the Viterbi forward/backtrack algorithm,
  the per-position coverage classification, and the SAM output format.
- The SAMv1 specification for the output record layout.
- Black-box behavior testing against `samtools targetcut 1.23.1`.

License: MIT OR Apache-2.0.  
Upstream credit: [samtools](https://github.com/samtools/samtools) (MIT/Expat).

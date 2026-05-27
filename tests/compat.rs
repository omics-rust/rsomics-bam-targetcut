//! Byte-exact compat against `samtools targetcut`.
//!
//! Runs both `samtools targetcut` and `rsomics-bam-targetcut` on two golden
//! fixtures and asserts identical stdout output line-by-line.
//!
//! Skipped gracefully if samtools is absent or older than 1.13.

use std::path::{Path, PathBuf};
use std::process::Command;

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bam-targetcut"))
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn samtools_version() -> Option<(u32, u32)> {
    let out = Command::new("samtools").arg("--version").output().ok()?;
    if !out.status.success() {
        return None;
    }
    let stdout = String::from_utf8_lossy(&out.stdout);
    let ver = stdout.lines().next()?.split_whitespace().nth(1)?;
    let mut it = ver.split('.');
    let major: u32 = it.next()?.parse().ok()?;
    let minor: u32 = it.next()?.parse().ok()?;
    Some((major, minor))
}

fn samtools_ready() -> bool {
    match samtools_version() {
        Some((maj, min)) if maj > 1 || (maj == 1 && min >= 13) => true,
        Some((maj, min)) => {
            eprintln!("SKIP targetcut compat: samtools {maj}.{min} (need >= 1.13)");
            false
        }
        None => {
            eprintln!("SKIP targetcut compat: samtools not found");
            false
        }
    }
}

fn run_samtools(bam: &Path, extra_args: &[&str]) -> Vec<u8> {
    let out = Command::new("samtools")
        .arg("targetcut")
        .args(extra_args)
        .arg(bam)
        .output()
        .expect("samtools targetcut failed");
    assert!(out.status.success(), "samtools targetcut exited non-zero");
    out.stdout
}

fn run_ours(bam: &Path, extra_args: &[&str]) -> Vec<u8> {
    let out = ours()
        .args(extra_args)
        .arg(bam)
        .output()
        .expect("rsomics-bam-targetcut failed");
    assert!(
        out.status.success(),
        "rsomics-bam-targetcut exited non-zero: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    out.stdout
}

fn compare(bam: &Path, samtools_args: &[&str], our_args: &[&str]) {
    let expected = run_samtools(bam, samtools_args);
    let actual = run_ours(bam, our_args);
    if expected != actual {
        // Show a readable diff on the first diverging line.
        let exp_lines: Vec<&str> = std::str::from_utf8(&expected)
            .unwrap_or("")
            .lines()
            .collect();
        let act_lines: Vec<&str> = std::str::from_utf8(&actual).unwrap_or("").lines().collect();
        for (i, (e, a)) in exp_lines.iter().zip(act_lines.iter()).enumerate() {
            if e != a {
                // Compare field-by-field (tab-split SAM).
                let ef: Vec<&str> = e.splitn(12, '\t').collect();
                let af: Vec<&str> = a.splitn(12, '\t').collect();
                for (fi, (ef, af)) in ef.iter().zip(af.iter()).enumerate() {
                    if ef != af {
                        panic!("Line {i} field {fi} differs\n  expected: {ef}\n  actual:   {af}");
                    }
                }
                panic!("Line {i} differs (possibly length):\n  expected: {e}\n  actual:   {a}");
            }
        }
        assert_eq!(
            exp_lines.len(),
            act_lines.len(),
            "line count differs: samtools {} vs ours {}",
            exp_lines.len(),
            act_lines.len()
        );
    }
}

#[test]
fn targetcut_default_large_fixture() {
    if !samtools_ready() {
        return;
    }
    let bam = golden("amplicon_large.bam");
    if !bam.exists() {
        eprintln!("SKIP: amplicon_large.bam not found");
        return;
    }
    compare(&bam, &[], &[]);
}

#[test]
fn targetcut_small_fixture_with_low_penalty() {
    if !samtools_ready() {
        return;
    }
    let bam = golden("amplicon.bam");
    if !bam.exists() {
        eprintln!("SKIP: amplicon.bam not found");
        return;
    }
    // Default penalty (14000) produces no output on this 200bp fixture;
    // use -i 100 to make both sides emit a target interval for comparison.
    compare(&bam, &["-i", "100"], &["-i", "100"]);
}

#[test]
fn targetcut_custom_emission_params() {
    if !samtools_ready() {
        return;
    }
    let bam = golden("amplicon_large.bam");
    if !bam.exists() {
        eprintln!("SKIP: amplicon_large.bam not found");
        return;
    }
    // Test that --em0/--em1/--em2 emission overrides propagate identically to
    // samtools' -0/-1/-2 flags. Negative emission values must be passed with =
    // to avoid clap misinterpreting them as flags.
    compare(
        &bam,
        &["-0", "-2", "-1", "2", "-2", "8"],
        &["--em0=-2", "--em1=2", "--em2=8"],
    );
}

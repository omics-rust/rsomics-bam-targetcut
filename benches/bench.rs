use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bam_targetcut(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bam-targetcut");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/amplicon.bam");
    c.bench_function("rsomics-bam-targetcut golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([bam.to_str().unwrap(), "-o", "/dev/null"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bam_targetcut);
criterion_main!(benches);

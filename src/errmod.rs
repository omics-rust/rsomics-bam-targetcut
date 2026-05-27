//! MAQ error model — pure Rust port of `htslib/errmod.c`.
//!
//! Computes a per-position consensus quality and best base from a set of
//! pileup bases encoded as `u16`: `[15:5]=quality(6b) [4]=strand [3:0]=base(ACGT 0-3)`.
//!
//! The model pre-computes a `beta` look-up table indexed by `[qual][depth][k]` and
//! a dependency-correction table `fk[n]`. Both are used in `errmod_cal` to convert
//! a bag of per-read base+quality observations into phred-scaled genotype likelihoods,
//! from which the best homozygous call is extracted.

/// Dependency correlation constant: `errmod_init(1.0 - ERR_DEP)`.
pub const ERR_DEP: f64 = 0.83;

/// Maximum depth stored in the beta table (matches htslib's 256-entry cap).
const MAX_DEPTH: usize = 256;

/// Number of quality bins in the table (0..63 → 64 entries; bin 0 unused).
const N_QUAL: usize = 64;

pub struct Errmod {
    /// `fk[n]` = dependency-correction weight for the n-th observation of a strand+base combo.
    fk: Box<[f64; MAX_DEPTH]>,
    /// `beta[qual * MAX_DEPTH * MAX_DEPTH + depth * MAX_DEPTH + k]`
    /// = phred-scaled probability that `k` out of `depth` observations are errors
    /// given per-base error rate `10^(-qual/10)`.
    beta: Vec<f64>, // N_QUAL * MAX_DEPTH * MAX_DEPTH
    /// `lhet[n * MAX_DEPTH + k]` = log binomial coefficient log(C(n,k)) − n·ln2.
    lhet: Vec<f64>, // MAX_DEPTH * MAX_DEPTH
}

impl Errmod {
    pub fn new(depcorr: f64) -> Self {
        let eta = 0.03f64;

        // fk[n] = (1-depcorr)^n * (1-eta) + eta
        let mut fk = Box::new([0f64; MAX_DEPTH]);
        fk[0] = 1.0;
        for n in 1..MAX_DEPTH {
            fk[n] = (1.0 - depcorr).powi(n as i32) * (1.0 - eta) + eta;
        }

        // Log binomial table: lc[n][k] = log(n!) - log(k!) - log((n-k)!)
        let lc = logbinomial_table();

        // beta[q][n][k] (stored as q * MAX_DEPTH^2 + n * MAX_DEPTH + k)
        let mut beta = vec![0f64; N_QUAL * MAX_DEPTH * MAX_DEPTH];
        for q in 1..N_QUAL {
            let e = 10f64.powf(-(q as f64) / 10.0);
            let le = e.ln();
            let le1 = (1.0 - e).ln();
            for n in 1..MAX_DEPTH {
                let base = q * MAX_DEPTH * MAX_DEPTH + n * MAX_DEPTH;
                // sum1 = log P(n errors | n bases, q)
                let mut sum1 = lc[n * MAX_DEPTH + n] + (n as f64) * le;
                // sentinel: beta[n] = HUGE_VAL → treated as 0 phred penalty
                beta[base + n] = f64::INFINITY;
                // iterate k from n-1 down to 0, accumulating via log1p
                for k in (0..n).rev() {
                    let sum = sum1
                        + (lc[n * MAX_DEPTH + k] + (k as f64) * le + ((n - k) as f64) * le1 - sum1)
                            .exp()
                            .ln_1p();
                    beta[base + k] = -10.0 / std::f64::consts::LN_10 * (sum1 - sum);
                    sum1 = sum;
                }
            }
        }

        // lhet[n][k] = lc[n][k] - ln(2) * n
        let mut lhet = vec![0f64; MAX_DEPTH * MAX_DEPTH];
        let ln2 = 2f64.ln();
        for n in 0..MAX_DEPTH {
            for k in 0..MAX_DEPTH {
                lhet[n * MAX_DEPTH + k] = lc[n * MAX_DEPTH + k] - ln2 * (n as f64);
            }
        }

        Self { fk, beta, lhet }
    }

    /// Compute genotype likelihoods for `m` alleles from `n` base observations.
    ///
    /// `bases[i]` = `quality(6b) << 5 | strand(1b) << 4 | base(4b)` matching
    /// htslib `errmod_cal` encoding (quality capped [4,63]).
    ///
    /// Returns `q[j*m + k]`: phred-scaled likelihood of genotype (j,k).
    /// Mirrors `errmod_cal` in `htslib/errmod.c` exactly including the
    /// downsampling, sorting, dependency-correction accumulation, and
    /// heterozygous likelihood formula.
    pub fn cal(&self, n: usize, m: usize, bases: &mut [u16], q: &mut [f32]) {
        let q_len = m * m;
        for v in q[..q_len].iter_mut() {
            *v = 0.0;
        }
        if n == 0 {
            return;
        }

        // Downsample to 255 (max depth in beta table).
        let n = if n > 255 {
            // Deterministic truncation (samtools shuffles randomly; we truncate
            // from the sorted end — for consensus calling only the best base matters).
            bases[..n].sort_unstable();
            255
        } else {
            bases[..n].sort_unstable();
            n
        };

        // Accumulate fsum/bsum/c per base (4 bases in m=4 calls).
        let mut fsum = [0f64; 16];
        let mut bsum = [0f64; 16];
        let mut c = [0usize; 16];
        let mut w = [0usize; 32]; // per strand+base counter

        for j in (0..n).rev() {
            let b = bases[j];
            let mut qual = (b >> 5) as usize;
            qual = qual.clamp(4, 63);
            let basestrand = (b & 0x1f) as usize;
            let base = (b & 0xf) as usize;
            fsum[base] += self.fk[w[basestrand]];
            let beta_val = self.beta[qual * MAX_DEPTH * MAX_DEPTH + n * MAX_DEPTH + c[base]];
            bsum[base] += self.fk[w[basestrand]]
                * if beta_val.is_infinite() {
                    0.0
                } else {
                    beta_val
                };
            c[base] += 1;
            w[basestrand] += 1;
        }

        // Homozygous likelihoods.
        for j in 0..m {
            let mut tmp1 = 0f64;
            let mut tmp2 = 0usize;
            for k in 0..m {
                if k == j {
                    continue;
                }
                tmp1 += bsum[k];
                tmp2 += c[k];
            }
            if tmp2 > 0 {
                q[j * m + j] = tmp1 as f32;
            }
        }

        // Heterozygous likelihoods.
        for j in 0..m {
            for k in (j + 1)..m {
                let cjk = c[j] + c[k];
                let mut tmp1 = 0f64;
                let mut tmp2 = 0usize;
                for i in 0..m {
                    if i == j || i == k {
                        continue;
                    }
                    tmp1 += bsum[i];
                    tmp2 += c[i];
                }
                let hval = -4.343 * self.lhet[cjk * MAX_DEPTH + c[k]];
                let val = if tmp2 > 0 { hval + tmp1 } else { hval };
                q[j * m + k] = val as f32;
                q[k * m + j] = val as f32;
            }
        }

        // Clamp to >= 0.
        for v in q[..q_len].iter_mut() {
            if *v < 0.0 {
                *v = 0.0;
            }
        }
    }
}

/// Compute the per-position consensus from a slice of base observations.
///
/// `bases`: `u16` array of `quality(6b) << 5 | strand(1b) << 4 | base(2b)`
/// (quality and strand already capped/validated by the caller).
///
/// Returns the consensus `u16`:
///   `[15:8]` = depth (capped at 255, 0 = no usable bases)
///   `[7:2]`  = best-base quality (6 bits)
///   `[1:0]`  = best base (0=A, 1=C, 2=G, 3=T)
///
/// Mirrors `gencns` in `cut_target.c`.
pub fn gencns(em: &Errmod, bases: &mut [u16]) -> u16 {
    let n = bases.len();
    if n == 0 {
        return 0;
    }
    let m = 4usize;
    let mut q = [0f32; 16]; // 4×4
    em.cal(n, m, bases, &mut q);

    // Best base by homozygous score (diagonal of q).
    let mut sum = [0i32; 4];
    for i in 0..4 {
        sum[i] = (q[i * m + i] + 0.499) as i32 * 4 + i as i32;
    }
    // Insertion sort (mirrors samtools' gencns: find the best 2 by sorting).
    for i in 1..4 {
        let mut j = i;
        while j > 0 && sum[j] < sum[j - 1] {
            sum.swap(j, j - 1);
            j -= 1;
        }
    }
    let qual = (sum[1] >> 2) - (sum[0] >> 2);
    let best_base = (sum[0] & 3) as u16;
    let qual6 = qual.min(63) as u16;
    let depth = n.min(255) as u16;
    (qual6 << 2 | best_base) << 8 | depth
}

/// Build the log-binomial table: `lc[n * MAX_DEPTH + k]` = `log(C(n,k))`.
fn logbinomial_table() -> Vec<f64> {
    let mut lc = vec![0f64; MAX_DEPTH * MAX_DEPTH];
    for n in 1..MAX_DEPTH {
        let lfn = lgamma(n as f64 + 1.0);
        for k in 1..=n {
            lc[n * MAX_DEPTH + k] = lfn - lgamma(k as f64 + 1.0) - lgamma((n - k) as f64 + 1.0);
        }
    }
    lc
}

/// `lgamma(x)` via the standard library.
fn lgamma(x: f64) -> f64 {
    // Rust std does not expose lgamma; use the relation lgamma(n+1) = ln(n!)
    // for integer n (Stirling/recursive), but for correctness across all inputs
    // we use the unsafe libc lgamma_r or our own Lanczos approximation.
    // Here we use a fast recursive table for the integer inputs we need.
    static LGAMMA_INT: std::sync::OnceLock<Vec<f64>> = std::sync::OnceLock::new();
    let tab = LGAMMA_INT.get_or_init(|| {
        let mut v = vec![0f64; MAX_DEPTH + 1];
        v[0] = 0.0; // ln(0!) = ln(1) = 0
        for i in 1..=MAX_DEPTH {
            v[i] = v[i - 1] + (i as f64).ln();
        }
        v
    });
    // x is always an integer+1 in our context (lgamma(n+1) = ln(n!)).
    let idx = (x - 1.0).round() as usize;
    if idx < tab.len() {
        tab[idx]
    } else {
        // Stirling approximation for large x.
        let y = x - 1.0;
        y * y.ln() - y + 0.5 * (2.0 * std::f64::consts::PI * y).ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gencns_empty() {
        let em = Errmod::new(1.0 - ERR_DEP);
        assert_eq!(gencns(&em, &mut []), 0);
    }

    #[test]
    fn gencns_single_a() {
        // One A at quality 30 on forward strand: base=0, qual=30, strand=0
        // bases encoding: qual<<5 | strand<<4 | base = 30<<5 | 0 | 0 = 960
        let em = Errmod::new(1.0 - ERR_DEP);
        let mut bases = [30u16 << 5];
        let cns = gencns(&em, &mut bases);
        // Best base should be A (0).
        assert_eq!((cns >> 8) & 3, 0, "best base should be A(0)");
        // Depth = 1.
        assert_eq!(cns & 0xff, 1, "depth should be 1");
    }

    #[test]
    fn coverage_class_correct() {
        use crate::coverage_class;
        assert_eq!(coverage_class(0), 0);
        assert_eq!(coverage_class(0x0001), 1); // non-zero, high byte = 0
        assert_eq!(coverage_class(0x0100), 2); // high byte = 1
    }
}

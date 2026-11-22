//! QC numbers on a count matrix: log-CPM, sample correlation and a two-component PCA.
//! Plain log-CPM values are used. They are not DESeq2 variance-stabilised values.
//! No file or network I/O.

use std::fmt;

use crate::counts::CountMatrix;

/// Why a QC number cannot be computed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QcError {
    /// The matrix has fewer samples than the analysis needs.
    TooFewSamples { need: usize, have: usize },
    /// The named sample has no counts in any gene.
    ZeroLibrary(String),
}

impl fmt::Display for QcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QcError::TooFewSamples { need, have } => {
                write!(f, "need at least {need} samples, have {have}")
            }
            QcError::ZeroLibrary(s) => {
                write!(f, "sample '{s}' has zero total counts and cannot be scaled")
            }
        }
    }
}

impl std::error::Error for QcError {}

fn need_samples(m: &CountMatrix, need: usize) -> Result<(), QcError> {
    let have = m.samples.len();
    if have < need {
        return Err(QcError::TooFewSamples { need, have });
    }
    Ok(())
}

/// `log2(count / library_size * 1e6 + 1)`, genes by samples. The library size of a sample is the
/// sum of its counts. A sample with library size 0 is `ZeroLibrary`.
pub fn log_cpm(m: &CountMatrix) -> Result<Vec<Vec<f64>>, QcError> {
    let n = m.samples.len();
    let mut lib = vec![0u64; n];
    for row in &m.counts {
        for (l, c) in lib.iter_mut().zip(row) {
            *l += c;
        }
    }
    if let Some(i) = lib.iter().position(|&l| l == 0) {
        return Err(QcError::ZeroLibrary(m.samples[i].clone()));
    }
    Ok(m.counts
        .iter()
        .map(|row| {
            row.iter()
                .zip(&lib)
                .map(|(&c, &l)| (c as f64 / l as f64 * 1e6 + 1.0).log2())
                .collect()
        })
        .collect())
}

/// Pearson correlation between samples on `log_cpm`, over genes with a total count above 0.
/// Needs 2 samples. The result is samples by samples. A sample with no variance over those genes
/// gives NaN against the others.
pub fn sample_correlation(m: &CountMatrix) -> Result<Vec<Vec<f64>>, QcError> {
    need_samples(m, 2)?;
    let l = log_cpm(m)?;
    let n = m.samples.len();
    let rows: Vec<&Vec<f64>> = l
        .iter()
        .zip(&m.counts)
        .filter(|(_, c)| c.iter().sum::<u64>() > 0)
        .map(|(r, _)| r)
        .collect();
    let g = rows.len() as f64;
    let mean: Vec<f64> = (0..n)
        .map(|j| rows.iter().map(|r| r[j]).sum::<f64>() / g)
        .collect();
    let mut cov = vec![vec![0.0; n]; n];
    for r in &rows {
        for i in 0..n {
            for j in i..n {
                cov[i][j] += (r[i] - mean[i]) * (r[j] - mean[j]);
            }
        }
    }
    let mut out = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i..n {
            let v = if i == j && cov[i][i] > 0.0 {
                1.0
            } else {
                (cov[i][j] / (cov[i][i] * cov[j][j]).sqrt()).clamp(-1.0, 1.0)
            };
            out[i][j] = v;
            out[j][i] = v;
        }
    }
    Ok(out)
}

/// The first two principal components.
#[derive(Debug, Clone, PartialEq)]
pub struct Pca {
    /// One `[PC1, PC2]` pair per sample, in matrix sample order.
    pub scores: Vec<[f64; 2]>,
    /// Share of the variance of the selected genes on PC1 and PC2 (0 to 1).
    pub variance_explained: [f64; 2],
}

/// PCA on `log_cpm` of the `top_genes` genes with the highest variance across samples (all genes
/// if the matrix has fewer). Each gene is centred over samples. Needs 3 samples.
pub fn pca(m: &CountMatrix, top_genes: usize) -> Result<Pca, QcError> {
    need_samples(m, 3)?;
    let l = log_cpm(m)?;
    let n = m.samples.len();
    let var = |r: &Vec<f64>| {
        let mean = r.iter().sum::<f64>() / n as f64;
        r.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / (n as f64 - 1.0)
    };
    let mut order: Vec<(usize, f64)> = l.iter().map(var).enumerate().collect();
    // Highest variance first; ties keep matrix order.
    order.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
    order.truncate(top_genes);
    // Gram matrix of the centred samples (samples by samples).
    let mut gram = vec![vec![0.0; n]; n];
    for &(gi, _) in &order {
        let r = &l[gi];
        let mean = r.iter().sum::<f64>() / n as f64;
        for i in 0..n {
            for j in i..n {
                gram[i][j] += (r[i] - mean) * (r[j] - mean);
            }
        }
    }
    for i in 0..n {
        for j in 0..i {
            gram[i][j] = gram[j][i];
        }
    }
    let total: f64 = (0..n).map(|i| gram[i][i]).sum();
    let (vals, vecs) = jacobi_eigen(gram);
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| vals[b].total_cmp(&vals[a]));
    let mut scores = vec![[0.0; 2]; n];
    let mut ve = [0.0; 2];
    for (k, &e) in idx.iter().take(2).enumerate() {
        let lambda = vals[e].max(0.0);
        // Fix the sign: the component with the largest magnitude is positive.
        let big = (0..n)
            .max_by(|&a, &b| vecs[a][e].abs().total_cmp(&vecs[b][e].abs()))
            .unwrap();
        let sign = if vecs[big][e] < 0.0 { -1.0 } else { 1.0 };
        for i in 0..n {
            scores[i][k] = sign * vecs[i][e] * lambda.sqrt();
        }
        ve[k] = if total > 0.0 { lambda / total } else { 0.0 };
    }
    Ok(Pca {
        scores,
        variance_explained: ve,
    })
}

/// Eigenvalues and eigenvectors (as columns) of a symmetric matrix, by cyclic Jacobi rotations.
fn jacobi_eigen(mut a: Vec<Vec<f64>>) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = a.len();
    let mut v = vec![vec![0.0; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| (0..n).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum();
        let diag: f64 = (0..n).map(|i| a[i][i] * a[i][i]).sum();
        if off <= 1e-30 * diag.max(1e-300) {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for row in v.iter_mut() {
                    let (vkp, vkq) = (row[p], row[q]);
                    row[p] = c * vkp - s * vkq;
                    row[q] = s * vkp + c * vkq;
                }
            }
        }
    }
    ((0..n).map(|i| a[i][i]).collect(), v)
}

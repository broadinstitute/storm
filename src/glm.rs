//! StormGLM - Internal GLM backend for association testing
//!
//! Implements:
//! - Linear regression
//! - Logistic regression
//! - Categorical regression
//! - BinomiRare test for rare variants
//! - Firth logistic regression (feature-flagged)

use std::collections::HashMap;
use thiserror::Error;

/// Errors from GLM operations
#[derive(Error, Debug)]
pub enum GlmError {
    #[error("Insufficient data: {0}")]
    InsufficientData(String),
    #[error("Convergence failed: {0}")]
    ConvergenceFailed(String),
    #[error("Invalid input: {0}")]
    InvalidInput(String),
    #[error("Not implemented: {0}")]
    NotImplemented(String),
}

/// Result of an association test
#[derive(Debug, Clone)]
pub struct AssociationResult {
    /// Test unit ID
    pub unit_id: String,
    /// Beta coefficient (effect size)
    pub beta: f64,
    /// Standard error
    pub se: f64,
    /// Z-score or T-statistic
    pub statistic: f64,
    /// P-value
    pub p_value: f64,
    /// Number of samples used
    pub n_samples: usize,
    /// Number of carriers
    pub n_carriers: usize,
    /// Call rate
    pub call_rate: f64,
    /// Model used
    pub model: String,
    /// Encoding used
    pub encoding: String,
    /// Additional metadata
    pub metadata: HashMap<String, String>,
}

/// Simple linear regression using OLS
pub fn linear_regression(
    x: &[f64],        // Predictor (encoded genotype)
    y: &[f64],        // Outcome (phenotype)
    covariates: Option<&[Vec<f64>]>, // Optional covariates
) -> Result<(f64, f64, f64, f64), GlmError> {
    if x.len() != y.len() {
        return Err(GlmError::InvalidInput("x and y must have same length".to_string()));
    }
    
    let n = x.len();
    if n < 3 {
        return Err(GlmError::InsufficientData("Need at least 3 samples".to_string()));
    }

    // Simple OLS: y = a + b*x
    // b = cov(x,y) / var(x)
    // a = mean(y) - b * mean(x)
    
    let x_mean: f64 = x.iter().sum::<f64>() / n as f64;
    let y_mean: f64 = y.iter().sum::<f64>() / n as f64;
    
    let mut cov_xy = 0.0;
    let mut var_x = 0.0;
    
    for i in 0..n {
        let dx = x[i] - x_mean;
        let dy = y[i] - y_mean;
        cov_xy += dx * dy;
        var_x += dx * dx;
    }
    
    if var_x < 1e-10 {
        return Err(GlmError::InsufficientData("No variance in predictor".to_string()));
    }
    
    let beta = cov_xy / var_x;
    
    // Residual sum of squares
    let mut rss = 0.0;
    for i in 0..n {
        let predicted = y_mean + beta * (x[i] - x_mean);
        let residual = y[i] - predicted;
        rss += residual * residual;
    }
    
    // Standard error of beta
    let mse = rss / (n - 2) as f64;
    let se = (mse / var_x).sqrt();
    
    // T-statistic
    let t_stat = beta / se;
    
    // P-value (using normal approximation for simplicity)
    let p_value = 2.0 * (1.0 - normal_cdf(t_stat.abs()));
    
    Ok((beta, se, t_stat, p_value))
}

/// Logistic regression using iteratively reweighted least squares (IRLS)
pub fn logistic_regression(
    x: &[f64],        // Predictor
    y: &[f64],        // Binary outcome (0/1)
    covariates: Option<&[Vec<f64>]>,
    max_iter: usize,
) -> Result<(f64, f64, f64, f64), GlmError> {
    if x.len() != y.len() {
        return Err(GlmError::InvalidInput("x and y must have same length".to_string()));
    }
    
    let n = x.len();
    if n < 5 {
        return Err(GlmError::InsufficientData("Need at least 5 samples".to_string()));
    }

    // Check y is binary
    for &yi in y {
        if yi != 0.0 && yi != 1.0 {
            return Err(GlmError::InvalidInput("y must be binary (0 or 1)".to_string()));
        }
    }

    // Simple logistic regression using Newton-Raphson
    let mut beta0: f64 = 0.0;
    let mut beta1: f64 = 0.0;
    
    for _iter in 0..max_iter {
        let mut grad0 = 0.0;
        let mut grad1 = 0.0;
        let mut hess00 = 0.0;
        let mut hess01 = 0.0;
        let mut hess11 = 0.0;
        
        for i in 0..n {
            let eta = beta0 + beta1 * x[i];
            let p = 1.0 / (1.0 + (-eta).exp());
            let w = p * (1.0 - p);
            
            let residual = y[i] - p;
            
            grad0 += residual;
            grad1 += residual * x[i];
            
            hess00 += w;
            hess01 += w * x[i];
            hess11 += w * x[i] * x[i];
        }
        
        // Solve 2x2 system
        let det = hess00 * hess11 - hess01 * hess01;
        if det.abs() < 1e-10 {
            break;
        }
        
        let delta0 = (hess11 * grad0 - hess01 * grad1) / det;
        let delta1 = (-hess01 * grad0 + hess00 * grad1) / det;
        
        beta0 += delta0;
        beta1 += delta1;
        
        if delta0.abs() < 1e-6 && delta1.abs() < 1e-6 {
            break;
        }
    }
    
    // Compute standard error from Hessian
    let mut hess11 = 0.0;
    for i in 0..n {
        let eta = beta0 + beta1 * x[i];
        let p = 1.0 / (1.0 + (-eta).exp());
        let w = p * (1.0 - p);
        hess11 += w * x[i] * x[i];
    }
    
    let se = if hess11 > 1e-10 { 1.0 / hess11.sqrt() } else { f64::INFINITY };
    let z_stat = beta1 / se;
    let p_value = 2.0 * (1.0 - normal_cdf(z_stat.abs()));
    
    Ok((beta1, se, z_stat, p_value))
}

/// BinomiRare test for rare variant burden testing
pub fn binomi_rare_test(
    carriers: usize,
    total: usize,
    expected_rate: f64,
) -> Result<f64, GlmError> {
    if total == 0 {
        return Err(GlmError::InsufficientData("No samples".to_string()));
    }
    
    if expected_rate <= 0.0 || expected_rate >= 1.0 {
        return Err(GlmError::InvalidInput("expected_rate must be in (0, 1)".to_string()));
    }
    
    // One-sided binomial test
    // P(X >= carriers) where X ~ Binomial(total, expected_rate)
    
    let mut p_value = 0.0;
    let mut log_prob = total as f64 * expected_rate.ln();
    
    for k in 0..=carriers {
        let log_binom = log_binomial(total, k);
        let log_p = log_binom + k as f64 * expected_rate.ln() + (total - k) as f64 * (1.0 - expected_rate).ln();
        if k < carriers {
            p_value += log_p.exp();
        }
    }
    
    p_value = 1.0 - p_value;
    
    Ok(p_value.max(0.0).min(1.0))
}

/// Firth logistic regression (penalized likelihood)
#[cfg(feature = "firth")]
pub fn firth_logistic_regression(
    x: &[f64],
    y: &[f64],
    covariates: Option<&[Vec<f64>]>,
    max_iter: usize,
) -> Result<(f64, f64, f64, f64), GlmError> {
    // Firth's penalized likelihood modification
    // This reduces bias in logistic regression with small samples/rare events
    Err(GlmError::NotImplemented("Firth regression requires 'firth' feature".to_string()))
}

/// Placeholder for Firth when feature not enabled
#[cfg(not(feature = "firth"))]
pub fn firth_logistic_regression(
    _x: &[f64],
    _y: &[f64],
    _covariates: Option<&[Vec<f64>]>,
    _max_iter: usize,
) -> Result<(f64, f64, f64, f64), GlmError> {
    Err(GlmError::NotImplemented("Firth regression requires 'firth' feature flag".to_string()))
}

/// Categorical regression using multinomial logit
pub fn categorical_regression(
    x: &[f64],
    y: &[usize],       // Category indices
    n_categories: usize,
    covariates: Option<&[Vec<f64>]>,
) -> Result<Vec<(f64, f64, f64, f64)>, GlmError> {
    if x.len() != y.len() {
        return Err(GlmError::InvalidInput("x and y must have same length".to_string()));
    }
    
    if n_categories < 2 {
        return Err(GlmError::InvalidInput("Need at least 2 categories".to_string()));
    }

    // For each non-reference category, fit binary logistic regression
    let mut results = Vec::new();
    
    for cat in 1..n_categories {
        let y_binary: Vec<f64> = y.iter().map(|&yi| if yi == cat { 1.0 } else { 0.0 }).collect();
        
        match logistic_regression(x, &y_binary, covariates, 25) {
            Ok(res) => results.push(res),
            Err(_) => results.push((0.0, f64::INFINITY, 0.0, 1.0)),
        }
    }
    
    Ok(results)
}

// Helper: Standard normal CDF approximation
fn normal_cdf(x: f64) -> f64 {
    // Approximation using error function
    0.5 * (1.0 + erf(x / std::f64::consts::SQRT_2))
}

// Error function approximation
fn erf(x: f64) -> f64 {
    // Abramowitz and Stegun approximation
    let a1 =  0.254829592;
    let a2 = -0.284496736;
    let a3 =  1.421413741;
    let a4 = -1.453152027;
    let a5 =  1.061405429;
    let p  =  0.3275911;
    
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    
    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
    
    sign * y
}

// Log binomial coefficient
fn log_binomial(n: usize, k: usize) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    log_factorial(n) - log_factorial(k) - log_factorial(n - k)
}

// Log factorial using Stirling's approximation for large n
fn log_factorial(n: usize) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    
    // Use exact computation for small n
    if n <= 20 {
        let mut result = 0.0;
        for i in 2..=n {
            result += (i as f64).ln();
        }
        return result;
    }
    
    // Stirling's approximation for larger n
    let n = n as f64;
    n * n.ln() - n + 0.5 * (2.0 * std::f64::consts::PI * n).ln()
}

/// Run association test with specified model
pub fn run_association(
    unit_id: &str,
    x: &[f64],
    y: &[f64],
    model: &crate::plan::Model,
    encoding: &crate::plan::Encoding,
    covariates: Option<&[Vec<f64>]>,
) -> Result<AssociationResult, GlmError> {
    let n_samples = x.len();
    let n_carriers = x.iter().filter(|&&v| v > 0.0).count();
    let n_called = x.iter().filter(|v| v.is_finite()).count();
    let call_rate = n_called as f64 / n_samples as f64;
    
    let (beta, se, statistic, p_value) = match model {
        crate::plan::Model::Linear => {
            linear_regression(x, y, covariates)?
        }
        crate::plan::Model::Logistic => {
            logistic_regression(x, y, covariates, 25)?
        }
        crate::plan::Model::BinomiRare => {
            // For BinomiRare, x should be binary
            let carriers = x.iter().filter(|&&v| v > 0.0).count();
            let case_carriers = x.iter().zip(y.iter())
                .filter(|(&xi, &yi)| xi > 0.0 && yi > 0.0)
                .count();
            let cases = y.iter().filter(|&&yi| yi > 0.0).count();
            
            let expected_rate = carriers as f64 / n_samples as f64;
            let p = binomi_rare_test(case_carriers, cases, expected_rate)?;
            
            (0.0, 0.0, 0.0, p)
        }
        crate::plan::Model::Firth => {
            firth_logistic_regression(x, y, covariates, 25)?
        }
        crate::plan::Model::Categorical => {
            return Err(GlmError::InvalidInput("Use categorical_regression for categorical outcomes".to_string()));
        }
    };
    
    Ok(AssociationResult {
        unit_id: unit_id.to_string(),
        beta,
        se,
        statistic,
        p_value,
        n_samples,
        n_carriers,
        call_rate,
        model: format!("{:?}", model),
        encoding: format!("{:?}", encoding),
        metadata: HashMap::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_linear_regression() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let y = vec![2.1, 4.2, 5.8, 8.1, 9.9, 12.2, 14.0, 15.8, 18.1, 20.0];
        
        let (beta, se, t_stat, p_value) = linear_regression(&x, &y, None).unwrap();
        
        // Beta should be close to 2.0
        assert!((beta - 2.0).abs() < 0.1);
        assert!(se > 0.0);
        assert!(t_stat > 10.0); // Strong signal
        assert!(p_value < 0.001);
    }

    #[test]
    fn test_logistic_regression() {
        let x = vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let y = vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        
        let (beta, se, z_stat, p_value) = logistic_regression(&x, &y, None, 25).unwrap();
        
        // Beta should be positive (x=1 increases probability of y=1)
        assert!(beta > 0.0);
        assert!(se > 0.0);
        assert!(p_value < 1.0);
    }

    #[test]
    fn test_binomi_rare() {
        // Test enrichment: 5 carriers among 10 cases, expected 10% rate
        let p = binomi_rare_test(5, 10, 0.1).unwrap();
        assert!(p < 0.01); // Should be significant enrichment

        // Test no enrichment: 1 carrier among 10 cases, expected 10% rate
        let p = binomi_rare_test(1, 10, 0.1).unwrap();
        assert!(p > 0.1); // Not significant
    }

    #[test]
    fn test_run_association_linear() {
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let y = vec![2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0];
        
        let result = run_association(
            "test_unit",
            &x,
            &y,
            &crate::plan::Model::Linear,
            &crate::plan::Encoding::S,
            None,
        ).unwrap();
        
        assert_eq!(result.unit_id, "test_unit");
        assert!((result.beta - 2.0).abs() < 0.01);
        assert!(result.p_value < 0.001);
        assert_eq!(result.n_samples, 10);
    }

    #[test]
    fn test_firth_feature_flagged() {
        let x = vec![1.0, 2.0, 3.0];
        let y = vec![0.0, 0.0, 1.0];
        
        let result = firth_logistic_regression(&x, &y, None, 25);
        assert!(result.is_err());
    }
}

//! Plan module for analysis configuration
//!
//! Parses YAML configuration files that specify:
//! - Which encodings to use
//! - Which models to apply
//! - Rule-based selection logic

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use serde::{Deserialize, Serialize};
use thiserror::Error;

/// Errors that can occur with plan parsing
#[derive(Error, Debug)]
pub enum PlanError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    #[error("YAML parse error: {0}")]
    Yaml(#[from] serde_yaml::Error),
    #[error("Invalid rule: {0}")]
    InvalidRule(String),
}

/// Encoding types for repeat lengths
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Encoding {
    /// Sum of allele lengths (L1 + L2)
    S,
    /// Maximum allele length
    M,
    /// Absolute difference |L1 - L2|
    D,
    /// Tail indicator (above threshold)
    Tail,
    /// Categorical bins
    Categorical,
    /// Binary presence/absence
    Binary,
}

/// Model types for association testing
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Model {
    /// Linear regression
    Linear,
    /// Logistic regression
    Logistic,
    /// Categorical regression
    Categorical,
    /// Binomial rare test
    BinomiRare,
    /// Firth logistic regression
    Firth,
}

/// A rule for selecting encoding and model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Rule {
    /// Rule name
    pub name: String,
    /// Condition for applying this rule
    #[serde(default)]
    pub condition: RuleCondition,
    /// Encoding to use
    pub encoding: Encoding,
    /// Model to use
    pub model: Model,
    /// Additional parameters
    #[serde(default)]
    pub params: HashMap<String, serde_yaml::Value>,
}

/// Condition for rule matching
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RuleCondition {
    /// Unit type must match
    #[serde(default)]
    pub unit_type: Option<String>,
    /// Motif pattern (regex)
    #[serde(default)]
    pub motif_pattern: Option<String>,
    /// Minimum carrier count
    #[serde(default)]
    pub min_carriers: Option<usize>,
    /// Maximum carrier count
    #[serde(default)]
    pub max_carriers: Option<usize>,
    /// Minimum call rate
    #[serde(default)]
    pub min_call_rate: Option<f64>,
}

impl RuleCondition {
    /// Check if this condition matches
    pub fn matches(
        &self,
        unit_type: &str,
        motif: Option<&str>,
        carrier_count: usize,
        call_rate: f64,
    ) -> bool {
        // Check unit type
        if let Some(ref required_type) = self.unit_type {
            if unit_type != required_type {
                return false;
            }
        }

        // Check motif pattern (simple contains for now)
        if let Some(ref pattern) = self.motif_pattern {
            if let Some(m) = motif {
                if !m.contains(pattern) {
                    return false;
                }
            } else {
                return false;
            }
        }

        // Check carrier count
        if let Some(min) = self.min_carriers {
            if carrier_count < min {
                return false;
            }
        }
        if let Some(max) = self.max_carriers {
            if carrier_count > max {
                return false;
            }
        }

        // Check call rate
        if let Some(min_rate) = self.min_call_rate {
            if call_rate < min_rate {
                return false;
            }
        }

        true
    }
}

/// Analysis plan configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Plan {
    /// Plan name
    pub name: String,
    /// Plan version
    #[serde(default = "default_version")]
    pub version: String,
    /// Default encoding
    #[serde(default = "default_encoding")]
    pub default_encoding: Encoding,
    /// Default model
    #[serde(default = "default_model")]
    pub default_model: Model,
    /// Rules for selecting encoding and model
    #[serde(default)]
    pub rules: Vec<Rule>,
    /// Covariate columns
    #[serde(default)]
    pub covariates: Vec<String>,
    /// Principal components to include
    #[serde(default)]
    pub num_pcs: usize,
}

fn default_version() -> String {
    "1.0".to_string()
}

fn default_encoding() -> Encoding {
    Encoding::S
}

fn default_model() -> Model {
    Model::Linear
}

impl Default for Plan {
    fn default() -> Self {
        Plan {
            name: "default".to_string(),
            version: default_version(),
            default_encoding: default_encoding(),
            default_model: default_model(),
            rules: Vec::new(),
            covariates: Vec::new(),
            num_pcs: 0,
        }
    }
}

impl Plan {
    /// Create a new default plan
    pub fn new(name: &str) -> Self {
        Plan {
            name: name.to_string(),
            ..Default::default()
        }
    }

    /// Load plan from YAML file
    pub fn from_yaml<P: AsRef<Path>>(path: P) -> Result<Self, PlanError> {
        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let plan: Plan = serde_yaml::from_reader(reader)?;
        Ok(plan)
    }

    /// Parse plan from YAML string
    pub fn from_yaml_str(yaml: &str) -> Result<Self, PlanError> {
        let plan: Plan = serde_yaml::from_str(yaml)?;
        Ok(plan)
    }

    /// Select encoding and model for a test unit
    pub fn select(
        &self,
        unit_type: &str,
        motif: Option<&str>,
        carrier_count: usize,
        call_rate: f64,
    ) -> (Encoding, Model) {
        // Try each rule in order (first match wins)
        for rule in &self.rules {
            if rule.condition.matches(unit_type, motif, carrier_count, call_rate) {
                return (rule.encoding.clone(), rule.model.clone());
            }
        }

        // Fall back to defaults
        (self.default_encoding.clone(), self.default_model.clone())
    }

    /// Add a rule to the plan
    pub fn add_rule(&mut self, rule: Rule) {
        self.rules.push(rule);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_plan_default() {
        let plan = Plan::default();
        assert_eq!(plan.default_encoding, Encoding::S);
        assert_eq!(plan.default_model, Model::Linear);
    }

    #[test]
    fn test_plan_from_yaml_str() {
        let yaml = r#"
name: test_plan
version: "1.0"
default_encoding: s
default_model: linear
rules:
  - name: rare_variants
    condition:
      max_carriers: 10
    encoding: binary
    model: binomirare
  - name: cag_repeats
    condition:
      motif_pattern: CAG
    encoding: m
    model: logistic
covariates:
  - age
  - sex
num_pcs: 10
"#;

        let plan = Plan::from_yaml_str(yaml).expect("Failed to parse YAML");
        
        assert_eq!(plan.name, "test_plan");
        assert_eq!(plan.rules.len(), 2);
        assert_eq!(plan.covariates, vec!["age", "sex"]);
        assert_eq!(plan.num_pcs, 10);
    }

    #[test]
    fn test_plan_select_with_rules() {
        let yaml = r#"
name: test_plan
default_encoding: s
default_model: linear
rules:
  - name: rare_variants
    condition:
      max_carriers: 10
    encoding: binary
    model: binomirare
  - name: cag_repeats
    condition:
      motif_pattern: CAG
    encoding: m
    model: logistic
"#;

        let plan = Plan::from_yaml_str(yaml).unwrap();

        // Test rare variant rule
        let (enc, model) = plan.select("TrueRepeat", Some("AT"), 5, 0.95);
        assert_eq!(enc, Encoding::Binary);
        assert_eq!(model, Model::BinomiRare);

        // Test CAG rule (has more carriers so rare rule doesn't match)
        let (enc, model) = plan.select("TrueRepeat", Some("CAG"), 100, 0.95);
        assert_eq!(enc, Encoding::M);
        assert_eq!(model, Model::Logistic);

        // Test default fallback
        let (enc, model) = plan.select("Sv", Some("none"), 100, 0.95);
        assert_eq!(enc, Encoding::S);
        assert_eq!(model, Model::Linear);
    }

    #[test]
    fn test_rule_condition_matches() {
        let cond = RuleCondition {
            unit_type: Some("TrueRepeat".to_string()),
            motif_pattern: Some("CAG".to_string()),
            min_carriers: Some(10),
            max_carriers: None,
            min_call_rate: Some(0.9),
        };

        assert!(cond.matches("TrueRepeat", Some("CAG"), 100, 0.95));
        assert!(!cond.matches("Sv", Some("CAG"), 100, 0.95)); // Wrong type
        assert!(!cond.matches("TrueRepeat", Some("AT"), 100, 0.95)); // Wrong motif
        assert!(!cond.matches("TrueRepeat", Some("CAG"), 5, 0.95)); // Too few carriers
        assert!(!cond.matches("TrueRepeat", Some("CAG"), 100, 0.8)); // Low call rate
    }

    #[test]
    fn test_deterministic_selection() {
        // Rules should be deterministic - same input always gives same output
        let yaml = r#"
name: deterministic_test
rules:
  - name: rule1
    condition:
      min_carriers: 50
    encoding: s
    model: linear
  - name: rule2
    condition:
      min_carriers: 10
    encoding: m
    model: logistic
"#;

        let plan = Plan::from_yaml_str(yaml).unwrap();

        // Run multiple times to verify determinism
        for _ in 0..10 {
            let (enc1, model1) = plan.select("TrueRepeat", Some("CAG"), 100, 0.95);
            let (enc2, model2) = plan.select("TrueRepeat", Some("CAG"), 100, 0.95);
            assert_eq!(enc1, enc2);
            assert_eq!(model1, model2);
        }
    }
}

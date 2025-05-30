#[cfg(feature = "python")]
use pyo3::prelude::*;

// Example library function
#[cfg(feature = "python")]
#[pyfunction]
pub fn add(a: i32, b: i32) -> i32 {
    a + b
}

/// Returns the version of the storm package as a string.
///
/// This function returns the version that was set in Cargo.toml at compile time.
#[cfg(feature = "python")]
#[pyfunction]
fn _version() -> PyResult<String> {
    Ok(env!("CARGO_PKG_VERSION").to_string())
}

#[cfg(feature = "python")]
#[pymodule]
fn storm(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(add, m)?)?;
    m.add_function(wrap_pyfunction!(_version, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[cfg(feature = "python")]
    fn test_add() {
        assert_eq!(add(2, 2), 4);
        assert_eq!(add(-1, 1), 0);
        assert_eq!(add(0, 0), 0);
    }

    #[test]
    #[cfg(feature = "python")]
    fn test_version() {
        let version = _version().unwrap();
        assert!(!version.is_empty());
        // Version should be in format x.y.z
        assert!(version.matches(r"^\d+\.\d+\.\d+$"));
    }
}

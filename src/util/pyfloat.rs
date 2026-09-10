//! Format `f64` values the way Python's `repr(float)` does, so that stats
//! files are byte-compatible with pairtools output.

use std::fmt::Write;

/// Format a float like Python's `repr`: shortest round-trip digits, fixed
/// notation for exponents in `[-4, 16)`, otherwise scientific notation with a
/// signed two-digit (minimum) exponent.
pub fn format_py_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x > 0.0 { "inf".into() } else { "-inf".into() };
    }
    if x == 0.0 {
        return if x.is_sign_negative() {
            "-0.0".into()
        } else {
            "0.0".into()
        };
    }
    // `{:e}` yields the shortest round-trip mantissa, e.g. "1.2345e-7".
    let sci = format!("{x:e}");
    let (mantissa, exp) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i32 = exp.parse().unwrap_or(0);
    let negative = mantissa.starts_with('-');
    let mantissa = mantissa.trim_start_matches('-');
    let digits: String = mantissa.chars().filter(|c| *c != '.').collect();
    let mut out = String::new();
    if negative {
        out.push('-');
    }
    if (-4..16).contains(&exp) {
        if exp >= 0 {
            let int_len = exp as usize + 1;
            if digits.len() <= int_len {
                out.push_str(&digits);
                for _ in digits.len()..int_len {
                    out.push('0');
                }
                out.push_str(".0");
            } else {
                out.push_str(&digits[..int_len]);
                out.push('.');
                out.push_str(&digits[int_len..]);
            }
        } else {
            out.push_str("0.");
            for _ in 0..(-exp - 1) {
                out.push('0');
            }
            out.push_str(&digits);
        }
    } else {
        out.push_str(&digits[..1]);
        if digits.len() > 1 {
            out.push('.');
            out.push_str(&digits[1..]);
        }
        let sign = if exp < 0 { '-' } else { '+' };
        let _ = write!(out, "e{sign}{:02}", exp.abs());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_python_repr() {
        assert_eq!(format_py_float(0.0), "0.0");
        assert_eq!(format_py_float(1.0), "1.0");
        assert_eq!(format_py_float(0.5), "0.5");
        assert_eq!(format_py_float(100.0), "100.0");
        assert_eq!(format_py_float(1234.5678), "1234.5678");
        assert_eq!(format_py_float(0.0001), "0.0001");
        assert_eq!(format_py_float(0.00001), "1e-05");
        assert_eq!(format_py_float(1e16), "1e+16");
        assert_eq!(format_py_float(1.5e16), "1.5e+16");
        assert_eq!(format_py_float(123456789012345.0), "123456789012345.0");
        assert_eq!(format_py_float(-2.5e-7), "-2.5e-07");
        assert_eq!(format_py_float(f64::NAN), "nan");
        assert_eq!(format_py_float(f64::INFINITY), "inf");
        assert_eq!(format_py_float(0.1 + 0.2), "0.30000000000000004");
        assert_eq!(format_py_float(1e22), "1e+22");
        assert_eq!(format_py_float(1e-100), "1e-100");
    }
}

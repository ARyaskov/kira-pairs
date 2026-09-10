//! Fast decimal integer parsing from byte slices.

/// Parse an unsigned decimal integer from ASCII bytes.
///
/// Accepts only `[0-9]+` (no sign, no whitespace). Returns `None` on empty
/// input, non-digit characters or overflow.
#[inline]
pub fn parse_u64(bytes: &[u8]) -> Option<u64> {
    if bytes.is_empty() || bytes.len() > 20 {
        return None;
    }
    let mut acc: u64 = 0;
    for &b in bytes {
        let d = b.wrapping_sub(b'0');
        if d > 9 {
            return None;
        }
        acc = acc.checked_mul(10)?.checked_add(u64::from(d))?;
    }
    Some(acc)
}

/// Parse a signed decimal integer (optional leading `-` or `+`).
#[inline]
pub fn parse_i64(bytes: &[u8]) -> Option<i64> {
    match bytes.first() {
        Some(b'-') => {
            let v = parse_u64(&bytes[1..])?;
            if v <= i64::MAX as u64 + 1 {
                Some((v as i64).wrapping_neg())
            } else {
                None
            }
        }
        Some(b'+') => i64::try_from(parse_u64(&bytes[1..])?).ok(),
        _ => i64::try_from(parse_u64(bytes)?).ok(),
    }
}

/// Parse a decimal float (`123`, `1.5`, `1e6`, `-2.5E-3`, `inf`, `nan`).
pub fn parse_f64(bytes: &[u8]) -> Option<f64> {
    let s = std::str::from_utf8(bytes).ok()?;
    s.parse::<f64>().ok()
}

/// Append the decimal representation of `v` to `out` without allocating.
#[inline]
pub fn write_u64(out: &mut Vec<u8>, mut v: u64) {
    let mut buf = [0u8; 20];
    let mut i = buf.len();
    loop {
        i -= 1;
        buf[i] = b'0' + (v % 10) as u8;
        v /= 10;
        if v == 0 {
            break;
        }
    }
    out.extend_from_slice(&buf[i..]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_unsigned() {
        assert_eq!(parse_u64(b"0"), Some(0));
        assert_eq!(parse_u64(b"123456789012"), Some(123_456_789_012));
        assert_eq!(parse_u64(b"18446744073709551615"), Some(u64::MAX));
        assert_eq!(parse_u64(b"18446744073709551616"), None);
        assert_eq!(parse_u64(b""), None);
        assert_eq!(parse_u64(b"12a"), None);
        assert_eq!(parse_u64(b"-1"), None);
        assert_eq!(parse_u64(b" 1"), None);
    }

    #[test]
    fn parses_signed() {
        assert_eq!(parse_i64(b"-5"), Some(-5));
        assert_eq!(parse_i64(b"+5"), Some(5));
        assert_eq!(parse_i64(b"9223372036854775807"), Some(i64::MAX));
        assert_eq!(parse_i64(b"-9223372036854775808"), Some(i64::MIN));
        assert_eq!(parse_i64(b"9223372036854775808"), None);
    }

    #[test]
    fn writes_unsigned() {
        let mut v = Vec::new();
        write_u64(&mut v, 0);
        v.push(b',');
        write_u64(&mut v, 1234567890123);
        assert_eq!(v, b"0,1234567890123");
    }
}

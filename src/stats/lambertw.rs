//! Principal branch of the Lambert W function on the real interval
//! `[-1/e, 0)`, as used by pairtools' naive library complexity estimate.

/// `W0(x)` for `-1/e <= x < 0` (returns NaN outside the supported range or
/// for positive inputs handled by the general iteration).
pub fn lambert_w0(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    let min = -(-1.0f64).exp();
    if x <= min {
        // Branch point (values slightly below -1/e are rounding noise).
        if x < min - 1e-12 {
            return f64::NAN;
        }
        return -1.0;
    }
    if x == 0.0 {
        return 0.0;
    }
    // Initial guess.
    let mut w = if x < -0.3 {
        // Series around the branch point.
        let p = (2.0 * (std::f64::consts::E * x + 1.0)).sqrt();
        -1.0 + p - p * p / 3.0 + 11.0 / 72.0 * p * p * p
    } else if x < 1.0 {
        x * (1.0 - x + 1.5 * x * x)
    } else {
        (1.0 + x).ln()
    };
    // Halley iteration.
    for _ in 0..64 {
        let ew = w.exp();
        let f = w * ew - x;
        let denom = ew * (w + 1.0) - (w + 2.0) * f / (2.0 * w + 2.0);
        if denom == 0.0 {
            break;
        }
        let delta = f / denom;
        w -= delta;
        if delta.abs() <= 1e-16 * (1.0 + w.abs()) {
            break;
        }
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_values() {
        assert!((lambert_w0(-0.1) - -0.11183255915896297).abs() < 1e-14);
        assert!((lambert_w0(-0.3) - -0.4894022271802149).abs() < 1e-14);
        assert!((lambert_w0(-0.36) - -0.8060843159708177).abs() < 1e-12);
        assert_eq!(lambert_w0(-(-1.0f64).exp()), -1.0);
        let w1 = lambert_w0(1.0);
        assert!((w1 - 0.5671432904097838).abs() < 1e-14, "{w1}");
        for i in 1..1000 {
            let x = -0.36787 * (i as f64) / 1000.0;
            let w = lambert_w0(x);
            assert!((w * w.exp() - x).abs() < 1e-13, "x={x} w={w}");
        }
    }
}

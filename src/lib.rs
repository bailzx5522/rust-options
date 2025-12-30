// lib.rs

pub mod pricing {
    use statrs::distribution::ContinuousCDF;
    use statrs::distribution::Normal;

    pub fn black_scholes_call_price(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let d1 = ((s / k).ln() + (r + 0.5 * v * v) * t) / (v * t.sqrt());
        let d2 = d1 - v * t.sqrt();
        let norm = Normal::new(0.0, 1.0).unwrap();
        s * norm.cdf(d1) - k * (-r * t).exp() * norm.cdf(d2)
    }

    pub fn black_scholes_put_price(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let d1 = ((s / k).ln() + (r + 0.5 * v * v) * t) / (v * t.sqrt());
        let d2 = d1 - v * t.sqrt();
        let norm = Normal::new(0.0, 1.0).unwrap();
        k * (-r * t).exp() * norm.cdf(-d2) - s * norm.cdf(-d1)
    }
}

pub mod greeks {
    use statrs::distribution::Continuous;
    use statrs::distribution::ContinuousCDF;
    use statrs::distribution::Normal;

    fn d1(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        ((s / k).ln() + (r + 0.5 * v * v) * t) / (v * t.sqrt())
    }

    fn d2(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        d1(s, k, t, r, v) - v * t.sqrt()
    }

    pub fn delta_call(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let d1 = ((s / k).ln() + (r + 0.5 * v * v) * t) / (v * t.sqrt());
        Normal::new(0.0, 1.0).unwrap().cdf(d1)
    }

    pub fn delta_put(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        delta_call(s, k, t, r, v) - 1.0
    }
    pub fn gamma(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let norm = Normal::new(0.0, 1.0).unwrap();
        norm.pdf(d1(s, k, t, r, v)) / (s * v * t.sqrt())
    }

    pub fn vega(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let norm = Normal::new(0.0, 1.0).unwrap();
        s * norm.pdf(d1(s, k, t, r, v)) * t.sqrt() / 100.0 // per 1% change
    }

    pub fn theta_call(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let norm = Normal::new(0.0, 1.0).unwrap();
        let d1 = d1(s, k, t, r, v);
        let d2 = d2(s, k, t, r, v);
        let term1 = -s * norm.pdf(d1) * v / (2.0 * t.sqrt());
        let term2 = r * k * (-r * t).exp() * norm.cdf(d2);
        (term1 - term2) / 365.0
    }

    pub fn theta_put(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let norm = Normal::new(0.0, 1.0).unwrap();
        let d1 = d1(s, k, t, r, v);
        let d2 = d2(s, k, t, r, v);
        let term1 = -s * norm.pdf(d1) * v / (2.0 * t.sqrt());
        let term2 = r * k * (-r * t).exp() * norm.cdf(-d2);
        (term1 + term2) / 365.0
    }

    pub fn rho_call(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let norm = Normal::new(0.0, 1.0).unwrap();
        k * t * (-r * t).exp() * norm.cdf(d2(s, k, t, r, v)) / 100.0
    }

    pub fn rho_put(s: f64, k: f64, t: f64, r: f64, v: f64) -> f64 {
        let norm = Normal::new(0.0, 1.0).unwrap();
        -k * t * (-r * t).exp() * norm.cdf(-d2(s, k, t, r, v)) / 100.0
    }

    // 使用期权价格反推IV
    pub fn implied_vol_call(
        binary_price: f64, // 市场报价，例如 0.62
        s: f64,            // 现价
        k: f64,            // 行权价
        t: f64,            // 到期时间（年）
    ) -> Option<f64> {
        let norm = Normal::new(0.0, 1.0).unwrap();

        // 第一步：反解出 d2
        let d2 = norm.inverse_cdf(binary_price);

        // 第二步：根据 d2 公式反解 sigma
        let ln_sk = (s / k).ln();
        let sqrt_t = t.sqrt();

        // d2 = [ln(S/K) + (r - q - σ²/2) T] / (σ √T)
        // 令 x = σ √T，则：
        // d2 * x = ln(S/K) + (r - q) T - 0.5 x²
        // 0.5 x² + d2 x - [ln(S/K) + (r - q)T] = 0

        let a = 0.5;
        let b = d2;
        let c = -(ln_sk);

        let discriminant = b * b - 4.0 * a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrt_disc = discriminant.sqrt();
        // 两个根，取正的且合理的那个（通常较小的那个）
        let x1 = (-b + sqrt_disc) / (2.0 * a);
        let x2 = (-b - sqrt_disc) / (2.0 * a);

        let x = if x1 > 0.0 { x1 } else { x2 };
        if x <= 0.0 {
            return None;
        }

        let sigma = x / sqrt_t;

        // 一般波动率不会超过 1000%，这里简单做个范围限制
        if sigma > 10.0 || sigma < 1e-6 {
            return None;
        }

        Some(sigma)
    }
}

pub mod strategy {
    /// Returns the payoff for a long call option
    pub fn long_call_payoff(s: f64, k: f64, premium: f64) -> f64 {
        (s - k).max(0.0) - premium
    }

    /// Returns the payoff for a short call option
    pub fn short_call_payoff(s: f64, k: f64, premium: f64) -> f64 {
        premium - (s - k).max(0.0)
    }

    /// Returns the payoff for a long put option
    pub fn long_put_payoff(s: f64, k: f64, premium: f64) -> f64 {
        (k - s).max(0.0) - premium
    }

    /// Returns the payoff for a short put option
    pub fn short_put_payoff(s: f64, k: f64, premium: f64) -> f64 {
        premium - (k - s).max(0.0)
    }

    /// Payoff for a long call spread (buy call at k1, sell call at k2)
    pub fn long_call_spread(s: f64, k1: f64, c1: f64, k2: f64, c2: f64) -> f64 {
        long_call_payoff(s, k1, c1) + short_call_payoff(s, k2, c2)
    }

    /// Payoff for a short call spread (sell call at k1, buy call at k2)
    pub fn short_call_spread(s: f64, k1: f64, c1: f64, k2: f64, c2: f64) -> f64 {
        short_call_payoff(s, k1, c1) + long_call_payoff(s, k2, c2)
    }

    /// Payoff for a long put spread (buy put at k1, sell put at k2)
    pub fn long_put_spread(s: f64, k1: f64, p1: f64, k2: f64, p2: f64) -> f64 {
        long_put_payoff(s, k1, p1) + short_put_payoff(s, k2, p2)
    }

    /// Payoff for a short put spread (sell put at k1, buy put at k2)
    pub fn short_put_spread(s: f64, k1: f64, p1: f64, k2: f64, p2: f64) -> f64 {
        short_put_payoff(s, k1, p1) + long_put_payoff(s, k2, p2)
    }

    /// Payoff for a long straddle (buy call and put at same strike)
    pub fn long_straddle(s: f64, k: f64, c: f64, p: f64) -> f64 {
        long_call_payoff(s, k, c) + long_put_payoff(s, k, p)
    }

    /// Payoff for a short straddle (sell call and put at same strike)
    pub fn short_straddle(s: f64, k: f64, c: f64, p: f64) -> f64 {
        short_call_payoff(s, k, c) + short_put_payoff(s, k, p)
    }

    /// Payoff for a long strangle (buy OTM call and OTM put)
    pub fn long_strangle(s: f64, kc: f64, c: f64, kp: f64, p: f64) -> f64 {
        long_call_payoff(s, kc, c) + long_put_payoff(s, kp, p)
    }

    /// Payoff for a short strangle (sell OTM call and OTM put)
    pub fn short_strangle(s: f64, kc: f64, c: f64, kp: f64, p: f64) -> f64 {
        short_call_payoff(s, kc, c) + short_put_payoff(s, kp, p)
    }
}

pub mod options {
    pub enum OptionType {
        Call,
        Put,
    }

    pub struct Option {
        pub option_type: OptionType,
        pub rfr: f64,
        pub strike: f64,
        pub spot: f64,
        pub iv: f64,
    }
}

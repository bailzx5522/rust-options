use statrs::distribution::Normal;
use statrs::distribution::ContinuousCDF;

/// Calculates the Black-Scholes price of a European call option
pub fn black_scholes_call_price(
    s: f64, // Spot price
    k: f64, // Strike price
    t: f64, // Time to maturity in years
    r: f64, // Risk-free rate
    v: f64  // Volatility
) -> f64 {
    let d1 = ((s / k).ln() + (r + 0.5 * v * v) * t) / (v * t.sqrt());
    let d2 = d1 - v * t.sqrt();
    let norm = Normal::new(0.0, 1.0).unwrap();
    s * norm.cdf(d1) - k * (-r * t).exp() * norm.cdf(d2)
}

pub fn simple_call_price(
    s: f64, // Spot price
    k: f64, // Strike price
    t: f64, // Time to maturity in years
    r: f64, // Risk-free rate
    v: f64  // Volatility
) -> f64 {
    let d1 = ((s / k).ln() + (r + 0.5 * v * v) * t) / (v * t.sqrt());
    let d2 = d1 - v * t.sqrt();
    let norm = Normal::new(0.0, 1.0).unwrap();
    let delta = Normal::new(0.0, 1.0).unwrap().cdf(d1);
    println!("delta:{}", delta);
    k * norm.cdf(d2)
}

/// Calculates the Black-Scholes price of a European put option
pub fn black_scholes_put_price(
    s: f64,
    k: f64,
    t: f64,
    r: f64,
    v: f64
) -> f64 {
    let d1 = ((s / k).ln() + (r + 0.5 * v * v) * t) / (v * t.sqrt());
    let d2 = d1 - v * t.sqrt();
    let norm = Normal::new(0.0, 1.0).unwrap();
    k * (-r * t).exp() * norm.cdf(-d2) - s * norm.cdf(-d1)
}

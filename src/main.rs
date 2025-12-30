mod options;
mod pricing;

use options::{Option, OptionType};
use pricing::black_scholes_call_price;

use crate::pricing::simple_call_price;


fn main() {

    let option = Option {
        option_type: OptionType::Call,
        rfr: 0.05,
        strike: 100.0,
        spot: 105.0,
        iv: 0.30
    };
    let t = 1.0/12.0;

    let option = Option {
        option_type: OptionType::Call,
        rfr: 0.0,
        strike: 2715.0,
        spot:   2720.0,
        // spot:   900.0,
        iv: 0.55
    };
    let t = 8.0 / (60 * 24 * 365 )as f64;

    // let price = black_scholes_call_price(option.spot, option.strike, t, option.rfr, option.iv);
    let price = simple_call_price(option.spot, option.strike, t, option.rfr, option.iv);
    println!("Call price {:.2}", price);
}

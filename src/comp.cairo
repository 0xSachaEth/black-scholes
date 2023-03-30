use option::OptionTrait;
use traits::Into;

use black_scholes::core::Fixed;
use black_scholes::core::FixedType;
use black_scholes::core::FixedImpl;
use black_scholes::core::FixedPartialOrd;


// PUBLIC

fn max (a: FixedType, b: FixedType) -> FixedType {
    if (a >= b) {
        return a;
    } else {
        return b;
    }
}

fn min (a: FixedType, b: FixedType) -> FixedType {
    if (a <= b) {
        return a;
    } else {
        return b;
    }
}

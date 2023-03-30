use black_scholes::core::Fixed;
use black_scholes::core::FixedType;
use black_scholes::core::ONE;
use black_scholes::core::ONE_u128;
use black_scholes::core::HALF;
use black_scholes::core::_felt_abs;
use black_scholes::core::_felt_sign;
use black_scholes::core::FixedInto;
use black_scholes::core::FixedPartialEq;
use black_scholes::core::FixedPartialOrd;
use black_scholes::core::FixedAdd;
use black_scholes::core::FixedAddEq;
use black_scholes::core::FixedSub;
use black_scholes::core::FixedSubEq;
use black_scholes::core::FixedMul;
use black_scholes::core::FixedMulEq;
use black_scholes::core::FixedDiv;
use black_scholes::trig::HALF_PI_u128;
use black_scholes::trig::PI_u128;
use traits::Into;
use gas::withdraw_gas;

fn exp(x : FixedType) -> FixedType {
    if(x.into() == 0){
        Fixed::from_unscaled_felt(1)
    } else {
        let mut acc = Fixed::from_unscaled_felt(1);  
        let mut n = 0;
        complete_taylor_series(acc, x, acc, n)
    }
}

fn complete_taylor_series(mut acc: FixedType, x: FixedType,mut last_value: FixedType, mut n: felt252) -> FixedType {
    loop {
 match gas::withdraw_gas_all(get_builtin_costs()) {
        Option::Some(_) => {
        },
        Option::None(_) => {
            let mut err_data = array::array_new();
            array::array_append(ref err_data, 'Out of gas');
            panic(err_data)
        },
    }
        if n == 51 {
            break acc;
        }
        n = n + 1;
        last_value = (last_value * x) / Fixed::from_unscaled_felt(n);  
        acc += last_value;
    }
}

#[test]
#[available_gas(9999999999999)]
fn test_exp() {
    let x_0 = Fixed::from_unscaled_felt(0);
    let res_0 = exp(x_0);
    assert(res_0 == Fixed::from_unscaled_felt(1), 'wtf 0');

    let x_1 = Fixed::from_unscaled_felt(1);
    let res_1 = exp(x_1);
    let l_b_1 = Fixed::from_felt(49990044073709551616); // = 2.71
    let h_b_1 = Fixed::from_felt(50290044073709551616); // = 2.72
    assert(res_1 > l_b_1, 'wtf l_b 1');
    assert(res_1 < h_b_1, 'wtf l_b 1');


    let x_2 = Fixed::from_unscaled_felt(2);
    let res_2 = exp(x_2);
    let l_b_2 = Fixed::from_felt(136135744073709551616); // = 7.38
    let h_b_2 = Fixed::from_felt(136605744073709551616); // = 7.40
    assert(res_2 > l_b_2, 'wtf l_b 2');
    assert(res_2 < h_b_2, 'wtf l_b 2');
}
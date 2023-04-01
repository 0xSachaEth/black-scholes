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
use debug::PrintTrait;


// Returns y, the exponent of x.
// Uses first 50 terms of taylor series expansion centered at 0.
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


// Used for below ln function approximation.
fn msb(x : FixedType) -> FixedType {
    match gas::withdraw_gas_all(get_builtin_costs()) {
        Option::Some(_) => {
        },
        Option::None(_) => {
            let mut err_data = array::array_new();
            array::array_append(ref err_data, 'Out of gas');
            panic(err_data)
        },
    }


    if(x <= Fixed::from_unscaled_felt(1)){
        Fixed::from_unscaled_felt(0)
    } else {
        let rest = msb(x / Fixed::from_unscaled_felt(2));
        rest + Fixed::from_unscaled_felt(1)
    }
}

// Returns y, the natural logarithm of x.
// Uses numerical approximation (Remez algorithm).
fn ln(x : FixedType) -> FixedType {
    match gas::withdraw_gas_all(get_builtin_costs()) {
        Option::Some(_) => {
        },
        Option::None(_) => {
            let mut err_data = array::array_new();
            array::array_append(ref err_data, 'Out of gas');
            panic(err_data)
        },
    }
    if(x.into() == ONE){
        Fixed::from_unscaled_felt(0)
    } else {    
        if(x < Fixed::from_unscaled_felt(1)){
             // ln(1/x) = -ln(x)
             ln( Fixed::from_unscaled_felt(1) / x );
             -x
        } else {
            let b = msb(x / Fixed::from_unscaled_felt(2));
            let norm = x / Fixed::from_unscaled_felt(2).pow(b);
            let d_1 = Fixed::from_felt(-1043548000000000000) * norm; // = -0.056570851
            let d_2 = (Fixed::from_felt(8249006700000000000) + d_1) * norm; // = 0.44717955
            let d_3 = (Fixed::from_felt(-27115917000000000000) + d_2) * norm; // = -1.4699568
            let d_4 = (Fixed::from_felt(52042002000000000000) + d_3) * norm; // = 2.8212026
            let d_5 = Fixed::from_felt(-32130426000000000000) + d_4; // = -1.7417939
            b * Fixed::from_felt(12786309000000000000) + d_5 // = ln(2) = 0.69314718
        }
    }
}



// Returns y, standard normal distribution at x.
// This computes e^(-x^2/2) / sqrt(2*pi).
fn std_normal(x : FixedType) -> FixedType{
    exp(-((x * x) / Fixed::from_unscaled_felt(2))) / Fixed::from_felt(46239130000000000000) // sqrt(2*pi) = 2.50662827463
}

// Returns y, cumulative normal distribution at x.
// Computed using a curve-fitting approximation.


fn std_normal_cdf(x : FixedType) -> FixedType {

    match gas::withdraw_gas_all(get_builtin_costs()) {
        Option::Some(_) => {
        },
        Option::None(_) => {
            let mut err_data = array::array_new();
            array::array_append(ref err_data, 'Out of gas');
            panic(err_data)
        },
    }

    if(x < Fixed::from_unscaled_felt(-5)){
        Fixed::from_unscaled_felt(0)
    } else {
        if(x > Fixed::from_unscaled_felt(5)){
            Fixed::from_unscaled_felt(1)
        } else {
            let b1 = Fixed::from_felt(5891549300000000000); // 0.31938153
            let b2 = Fixed::from_felt(-6577440800000000000); //  -0.356563782
            let b3 = Fixed::from_felt(32862468000000000000); // 1.781477937
            let b4 = Fixed::from_felt(-33596243000000000000); // -1.821255978
            let b5 = Fixed::from_felt(24539232000000000000); //  1.330274429
            let p = Fixed::from_felt(4273038800000000000); // 0.2316419;
            let c2 = Fixed::from_felt(7359186500000000000); // 0.3989423

            let abs_x = x;
            if(x < Fixed::from_unscaled_felt(0)){
                let abs_x = -x;
            }

            let t = Fixed::from_unscaled_felt(1) / (Fixed::from_unscaled_felt(1) + abs_x * p);
            let b = c2 / exp((x*x) / Fixed::from_unscaled_felt(2));
            let n = Fixed::from_unscaled_felt(1) - b*((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
            
            if(x < Fixed::from_unscaled_felt(0)){
                Fixed::from_unscaled_felt(1) - n
            } else {
                n
            }
        }
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
    assert(res_1 < h_b_1, 'wtf h_b 1');


    let x_2 = Fixed::from_unscaled_felt(2);
    let res_2 = exp(x_2);
    let l_b_2 = Fixed::from_felt(136135744073709551616); // = 7.38
    let h_b_2 = Fixed::from_felt(136605744073709551616); // = 7.40
    assert(res_2 > l_b_2, 'wtf l_b 2');
    assert(res_2 < h_b_2, 'wtf h_b 2');
}

#[test]
#[available_gas(9999999999999)]
fn test_ln() {
    let x_0 = Fixed::from_unscaled_felt(1);
    let res_0 = ln(x_0);
    assert(res_0 == Fixed::from_unscaled_felt(0), 'wtf 0');

    let x_1 = exp(Fixed::from_unscaled_felt(1));
    let res_1 = ln(x_1);
    let l_b_1 = Fixed::from_felt(18262277000000000000); // = 0.99
    let h_b_1 = Fixed::from_felt(18631212000000000000); // = 1.01
    assert(res_1 > l_b_1, 'wtf l_b 1');
    assert(res_1 < h_b_1, 'wtf h_b 1');


    let x_2 = exp(Fixed::from_unscaled_felt(2));
    let res_2 = ln(x_2);
    let l_b_2 = Fixed::from_felt(36709021000000000000); // = 1.99
    let h_b_2 = Fixed::from_felt(37077956000000000000); // = 2.01
    assert(res_2 > l_b_2, 'wtf l_b 2');
    assert(res_2 < h_b_2, 'wtf h_b 2');
}

#[test]
#[available_gas(9999999999999)]
fn test_std_normal() {
    let x_1 = Fixed::from_unscaled_felt(0);
    let res_1 = std_normal(x_1);
    let l_b_1 = Fixed::from_felt(7341804100000000000); // = 0.398
    let h_b_1 = Fixed::from_felt(7360250900000000000); // = 0.399
    assert(res_1 > l_b_1, 'wtf l_b 1');
    assert(res_1 < h_b_1, 'wtf h_b 1');

    let x_2 = Fixed::from_unscaled_felt(2);
    let res_2 = std_normal(x_2);
    let l_b_2 = Fixed::from_felt(977677440000000000); // = 0.053
    let h_b_2 = Fixed::from_felt(996124180000000000); // = 0.054
    assert(res_2 > l_b_2, 'wtf l_b 2');
    assert(res_2 < h_b_2, 'wtf h_b 2');
}

#[test]
#[available_gas(9999999999999)]
fn test_std_normal_cdf() {
    let x_1 = Fixed::from_unscaled_felt(-6);
    let res_1 = std_normal_cdf(x_1);
    assert(res_1 == Fixed::from_unscaled_felt(0), 'NOT 0');

    let x_2 = Fixed::from_unscaled_felt(6);
    let res_2 = std_normal_cdf(x_2);
    assert(res_2 == Fixed::from_unscaled_felt(1), 'NOT 1');

    let x_3 = Fixed::from_unscaled_felt(0);
    let res_3 = std_normal_cdf(x_3);
    res_3.into().print();

    // let x_4 = Fixed::from_felt(46116860000000000000);
    // let res_4 = std_normal_cdf(x_4);
    // res_4.into().print();

    let x_5 = Fixed::from_unscaled_felt(-1);
    let res_5 = std_normal_cdf(x_5);
    res_5.into().print();
}

// 15520070899922655376
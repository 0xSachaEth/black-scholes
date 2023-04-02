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

            let abs_x = absolute(x);

            let t = Fixed::from_unscaled_felt(1) / (Fixed::from_unscaled_felt(1) + abs_x * p);
            let b = c2 / exp((abs_x*abs_x) / Fixed::from_unscaled_felt(2));
            let n = Fixed::from_unscaled_felt(1) - b*((((b5*t+b4)*t+b3)*t+b2)*t+b1)*t;
            
            if(x < Fixed::from_unscaled_felt(0)){
                Fixed::from_unscaled_felt(1) - n
            } else {
                n
            }
        }
    }
}

fn absolute(x : FixedType) -> FixedType {
    if(x < Fixed::from_unscaled_felt(0)){
        - x
    } else {
        x
    }
} 

// Returns the internal Black-Scholes coefficients.
fn d1d2(t_annualised : FixedType, 
        volatility : FixedType, 
        spot : FixedType, 
        strike : FixedType, 
        rate : FixedType) -> (FixedType, FixedType) {
    let d1 = (ln(spot / strike) + (((volatility * volatility) / Fixed::from_unscaled_felt(2) + rate) * t_annualised)) / (volatility * t_annualised.sqrt());
    (d1, d1 - (volatility * t_annualised.sqrt()))
} 


// Returns the option's call and put delta value.
fn delta(t_annualised : FixedType, 
        volatility : FixedType, 
        spot : FixedType, 
        strike : FixedType, 
        rate : FixedType) -> (FixedType, FixedType) {
    let (d1, d2) = d1d2(t_annualised, volatility, spot, strike, rate);
    let call_delta = std_normal_cdf(d1);
    (call_delta, call_delta - Fixed::from_unscaled_felt(1))
} 

// Returns the option's gamma value.
fn gamma(t_annualised : FixedType, 
        volatility : FixedType, 
        spot : FixedType, 
        strike : FixedType, 
        rate : FixedType) -> FixedType {
    let (d1, d2) = d1d2(t_annualised, volatility, spot, strike, rate);
    std_normal(d1) / (spot * (volatility * t_annualised.sqrt()))
} 

//Returns the option's vega value.
fn vega(t_annualised : FixedType, 
        volatility : FixedType, 
        spot : FixedType, 
        strike : FixedType, 
        rate : FixedType) -> FixedType {
    let (d1, d2) = d1d2(t_annualised, volatility, spot, strike, rate);
    (t_annualised.sqrt() * (std_normal(d1) * spot)) / Fixed::from_unscaled_felt(100)
} 

// Returns the option's call and put rho value.
fn rho(t_annualised : FixedType, 
        volatility : FixedType, 
        spot : FixedType, 
        strike : FixedType, 
        rate : FixedType) -> (FixedType, FixedType) {
    let (d1, d2) = d1d2(t_annualised, volatility, spot, strike, rate);
    let strike_t = strike * t_annualised;
    let r_t = rate * t_annualised;
    let exp_trm = exp(- r_t);
    let lhs = strike_t * exp_trm;
    let d2_cdf = std_normal_cdf(d2);
    let d2_cdf_neg = std_normal_cdf(-d2);
    (lhs * d2_cdf / Fixed::from_unscaled_felt(100), lhs * d2_cdf_neg / Fixed::from_unscaled_felt(-100))
} 

// Returns the option's call and put theta value
fn theta(t_annualised : FixedType, 
        volatility : FixedType, 
        spot : FixedType, 
        strike : FixedType, 
        rate : FixedType) -> (FixedType, FixedType) {
    let (d1, d2) = d1d2(t_annualised, volatility, spot, strike, rate);
    let c1 = strike * rate;
    let c2 = c1 * exp(-(rate * t_annualised));
    let c3 = c2 * std_normal_cdf(d2);
    let p3 = c2 * std_normal_cdf(-d2);
    let c4 = ((spot * volatility) / (Fixed::from_unscaled_felt(2) * t_annualised.sqrt()));
    let c5 = std_normal(d1) * c4;
    ((-c5 -c3) / Fixed::from_unscaled_felt(365), (-c5 + p3) / Fixed::from_unscaled_felt(365))
} 


// Returns the option's call and put theta value
fn option_prices(t_annualised : FixedType, 
        volatility : FixedType, 
        spot : FixedType, 
        strike : FixedType, 
        rate : FixedType) -> (FixedType, FixedType) {
    let (d1, d2) = d1d2(t_annualised, volatility, spot, strike, rate);
    let exponent_term = exp(-rate * t_annualised);
    let strike_pv = exponent_term * strike;
    let spot_nd1 = spot * std_normal_cdf(d1);
    let strike_nd2 = strike_pv * std_normal_cdf(d2);
    (spot_nd1 - strike_nd2, (spot_nd1 - strike_nd2) + strike_pv - spot)
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
    let l_b_3 = Fixed::from_felt(9038904500000000000); // = 0.49
    let h_b_3 = Fixed::from_felt(9407839890450000000); // = 0.51
    assert(res_3 > l_b_3, 'wtf l_b 3');
    assert(res_3 < h_b_3, 'wtf h_b 3');

    let x_4 = Fixed::from_unscaled_felt(-1);
    let res_4 = std_normal_cdf(x_4);
    let l_b_4 = Fixed::from_felt(2767073173786896240); // = 0.15
    let h_b_4 = Fixed::from_felt(2951473173786896240); // = 0.16
    assert(res_4 > l_b_4, 'wtf l_b 4');
    assert(res_4 < h_b_4, 'wtf h_b 4');

    let x_5 = Fixed::from_unscaled_felt(1);
    let res_5 = std_normal_cdf(x_5);
    let l_b_5 = Fixed::from_felt(15495000000000000000); // = 0.84
    let h_b_5 = Fixed::from_felt(15679000000000000000); // = 0.85
    assert(res_5 > l_b_5, 'wtf l_b 4');
    assert(res_5 < h_b_5, 'wtf h_b 4');
}

#[test]
#[available_gas(9999999999999)]
fn test_d1d2() {
    let t_annualised = Fixed::from_unscaled_felt(1); // 1 year 
    let volatility = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(2); // 50%
    let spot = Fixed::from_unscaled_felt(50); // 
    let strike = Fixed::from_unscaled_felt(30);
    let rate = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(20); // 5%
    let (d_1, d_2) = d1d2(t_annualised, volatility, spot, strike, rate);
    let l_b_1 = Fixed::from_felt(25272000000000000000); // = 1.37
    let h_b_1 = Fixed::from_felt(25456500000000000000); // = 1.38
    let l_b_2 = Fixed::from_felt(16048660000000000000); // = 0.87   
    let h_b_2 = Fixed::from_felt(16233000000000000000); // = 0.88
    assert(d_1 > l_b_1, 'wtf l_b 1');
    assert(d_1 < h_b_1, 'wtf h_b 1');
    assert(d_2 > l_b_2, 'wtf l_b 2');
    assert(d_2 < h_b_2, 'wtf h_b 2');
}

#[test]
#[available_gas(9999999999999)]
fn test_delta() {
    let t_annualised = Fixed::from_unscaled_felt(1); // 1 year 
    let volatility = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(2); // 50%
    let spot = Fixed::from_unscaled_felt(50); // 
    let strike = Fixed::from_unscaled_felt(30);
    let rate = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(20); // 5%
    let (call_delta, put_delta) = delta(t_annualised, volatility, spot, strike, rate);
    let l_b_1 = Fixed::from_felt(16786000000000000000); // = 0.91
    let h_b_1 = Fixed::from_felt(16971000000000000000); // = 0.92
    let l_b_2 = Fixed::from_felt(-1586400000000000000); // = -0.086   
    let h_b_2 = Fixed::from_felt(-1567900000000000000); // = -0.085
    assert(call_delta > l_b_1, 'wtf l_b 1');
    assert(call_delta < h_b_1, 'wtf h_b 1');
    assert(put_delta > l_b_2, 'wtf l_b 2');
    assert(put_delta < h_b_2, 'wtf h_b 2');
}

#[test]
#[available_gas(9999999999999)]
fn test_gamma() {
    let t_annualised = Fixed::from_unscaled_felt(1); // 1 year 
    let volatility = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(2); // 50%
    let spot = Fixed::from_unscaled_felt(50); 
    let strike = Fixed::from_unscaled_felt(30);
    let rate = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(20); // 5%
    let gamma = gamma(t_annualised, volatility, spot, strike, rate);
    let l_b_1 = Fixed::from_felt(110680000000000000); // = 0.0060
    let h_b_1 = Fixed::from_felt(116214000000000000); // = 0.0063
    assert(gamma > l_b_1, 'wtf l_b 1');
    assert(gamma < h_b_1, 'wtf h_b 1');
}

#[test]
#[available_gas(9999999999999)]
fn test_vega() {
    let t_annualised = Fixed::from_unscaled_felt(1); // 1 year 
    let volatility = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(2); // 50%
    let spot = Fixed::from_unscaled_felt(50); 
    let strike = Fixed::from_unscaled_felt(30);
    let rate = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(20); // 5%
    let vega = vega(t_annualised, volatility, spot, strike, rate);
    let l_b_1 = Fixed::from_felt(1420390000000000000); // = 0.077
    let h_b_1 = Fixed::from_felt(1438890000000000000); // = 0.078
    assert(vega > l_b_1, 'wtf l_b 1');
    assert(vega < h_b_1, 'wtf h_b 1');
    // vega.into().print()
}

#[test]
#[available_gas(9999999999999)]
fn test_rho() {
    let t_annualised = Fixed::from_unscaled_felt(1); // 1 year 
    let volatility = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(2); // 50%
    let spot = Fixed::from_unscaled_felt(50); 
    let strike = Fixed::from_unscaled_felt(30);
    let rate = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(20); // 5%
    let (call_rho, put_rho) = rho(t_annualised, volatility, spot, strike, rate);
    let l_b_1 = Fixed::from_felt(4242750000000000000); // = 0.230
    let h_b_1 = Fixed::from_felt(4261150000000000000); // = 0.231
    let l_b_2 = Fixed::from_felt(-1014570000000000000); // = -0.055
    let h_b_2 = Fixed::from_felt(-996124000000000000); // = -0.054
    assert(call_rho > l_b_1, 'wtf l_b 1');
    assert(call_rho < h_b_1, 'wtf h_b 1');
    assert(put_rho > l_b_2, 'wtf l_b 2');
    assert(put_rho < h_b_2, 'wtf h_b 2');
}

#[test]
#[available_gas(9999999999999)]
fn test_theta() {
    let t_annualised = Fixed::from_unscaled_felt(1); // 1 year 
    let volatility = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(2); // 50%
    let spot = Fixed::from_unscaled_felt(50); 
    let strike = Fixed::from_unscaled_felt(30);
    let rate = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(20); // 5%
    let (call_theta, put_theta) = theta(t_annualised, volatility, spot, strike, rate);

    call_theta.into().print();
    let l_b_1 = Fixed::from_felt(-156797300000000000); // = -0.0085
    let h_b_1 = Fixed::from_felt(-154952600000000000); // = -0.0084

    let l_b_2 = Fixed::from_felt(-84855300000000000); // = -0.0046
    let h_b_2 = Fixed::from_felt(-83010300000000000); // = -0.0045
    assert(call_theta > l_b_1, 'wtf l_b 1');
    assert(call_theta < h_b_1, 'wtf h_b 1');
    assert(put_theta > l_b_2, 'wtf l_b 2');
    assert(put_theta < h_b_2, 'wtf h_b 2');
}

#[test]
#[available_gas(9999999999999)]
fn test_option_prices() {
    let t_annualised = Fixed::from_unscaled_felt(1); // 1 year 
    let volatility = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(2); // 50%
    let spot = Fixed::from_unscaled_felt(50); 
    let strike = Fixed::from_unscaled_felt(30);
    let rate = Fixed::from_unscaled_felt(1) / Fixed::from_unscaled_felt(20); // 5%
    let (call_price, put_price) = option_prices(t_annualised, volatility, spot, strike, rate);

    let l_b_1 = Fixed::from_felt(416890000000000000000); // = 22.6
    let h_b_1 = Fixed::from_felt(418741000000000000000); // = 22.7
    assert(call_price > l_b_1, 'wtf l_b 1');
    assert(call_price < h_b_1, 'wtf h_b 1');

    let l_b_2 = Fixed::from_felt(22320560000000000000); // = 1.21
    let h_b_2 = Fixed::from_felt(22505000000000000000); // = 1.22
    assert(put_price > l_b_2, 'wtf l_b 1');
    assert(put_price < h_b_2, 'wtf h_b 1');

}
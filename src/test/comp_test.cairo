// use option::OptionTrait;
// use traits::Into;

// use cubit::core::Fixed;
// use cubit::core::FixedInto;
// use cubit::core::FixedPartialEq;
// use cubit::comp::max;
// use cubit::comp::min;


// #[test]
// fn test_max() {
//     let a = Fixed::from_unscaled_felt(1);
//     let b = Fixed::from_unscaled_felt(0);
//     let c = Fixed::from_unscaled_felt(-1);

//     assert(max(a, a) == a, 'max(a, a)');
//     assert(max(a, b) == a, 'max(a, b)');
//     assert(max(a, c) == a, 'max(a, c)');

//     assert(max(b, a) == a, 'max(b, a)');
//     assert(max(b, b) == b, 'max(b, b)');
//     assert(max(b, c) == b, 'max(b, c)');

//     assert(max(c, a) == a, 'max(c, a)');
//     assert(max(c, b) == b, 'max(c, b)');
//     assert(max(c, c) == c, 'max(c, c)');
// }

// #[test]
// fn test_min() {
//     let a = Fixed::from_unscaled_felt(1);
//     let b = Fixed::from_unscaled_felt(0);
//     let c = Fixed::from_unscaled_felt(-1);

//     assert(min(a, a) == a, 'min(a, a)');
//     assert(min(a, b) == b, 'min(a, b)');
//     assert(min(a, c) == c, 'min(a, c)');

//     assert(min(b, a) == b, 'min(b, a)');
//     assert(min(b, b) == b, 'min(b, b)');
//     assert(min(b, c) == c, 'min(b, c)');

//     assert(min(c, a) == c, 'min(c, a)');
//     assert(min(c, b) == c, 'min(c, b)');
//     assert(min(c, c) == c, 'min(c, c)');
// }

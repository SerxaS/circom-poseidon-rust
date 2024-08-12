// Poseidon Hash implementation from Circom Library
// We have constants till t = 3

use halo2::{arithmetic::Field, halo2curves::bn256::Fr};

// Exp of S-box.
pub fn sigma(item: Fr) -> Fr {
    let item_2 = item * item;
    let item_4 = item_2 * item_2;
    item * item_4
}

// Adds round constants.
pub fn ark(t: usize, c: Vec<&'static str>, r: usize) {}

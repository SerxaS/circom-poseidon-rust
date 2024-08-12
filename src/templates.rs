// Poseidon Hash implementation from Circom Library
// We use parameters just for t = 3! Check here for more constants
// https://github.com/iden3/circomlib/blob/master/circuits/poseidon_constants.circom

use halo2::halo2curves::bn256::Fr;

/// Exp of S-box.
pub fn sigma(item: Fr) -> Fr {
    let item_2 = item * item;
    let item_4 = item_2 * item_2;
    item * item_4
}

/// Adds round constants.
pub fn ark(inputs: [Fr; 3], t: usize, c: Vec<Fr>, r: usize) -> [Fr; 3] {
    let mut state = [Fr::zero(); 3];

    for i in 0..t {
        state[i] = inputs[i] + c[i + r]
    }

    state
}

pub fn mix(state: [Fr; 3], t: usize, m: [[Fr; 3]; 3]) -> [Fr; 3] {
    let mut new_state = [Fr::zero(); 3];

    for i in 0..t {
        for j in 0..t {
            new_state[i] += m[j][i] * state[j];
        }
    }

    new_state
}

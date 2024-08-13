// Poseidon Hash implementation from Circom Library
// https://github.com/iden3/circomlib/blob/master/circuits/poseidon.circom
// We use parameters just for ninput = 2(t = 3)! Check here for more constants
// https://github.com/iden3/circomlib/blob/master/circuits/poseidon_constants.circom

use halo2::halo2curves::bn256::Fr;

use crate::utils::{poseidon_c, poseidon_m, poseidon_p, poseidon_s};

/// Exp of S-box.
pub fn sigma(item: Fr) -> Fr {
    let item_2 = item * item;
    let item_4 = item_2 * item_2;
    item * item_4
}

/// Adds round constants.
pub fn ark(inputs: [Fr; 2], t: usize, c: Vec<Fr>, r: usize) -> Vec<Fr> {
    let mut state = vec![Fr::zero(); t];
    for i in 0..t {
        state[i] = inputs[i] + c[i + r]
    }

    state
}

/// Computes MDS matrix for MixLayer operation.
pub fn mix(state: [Fr; 3], t: usize, m: [[Fr; 3]; 3]) -> [Fr; 3] {
    let mut new_state = [Fr::zero(); 3];

    for i in 0..t {
        for j in 0..t {
            new_state[i] += m[j][i] * state[j];
        }
    }

    new_state
}

/// Hash arithmetics.
pub fn poseidon_ex(inputs: [Fr; 2], initial_state: Fr, n_outs: usize) -> Fr {
    let n_rounds_p = [
        56, 57, 56, 60, 60, 63, 64, 63, 60, 66, 60, 65, 70, 60, 64, 68,
    ];
    let t = inputs.len() + 1;
    let n_rounds_f = 8;
    let n_round_p = n_rounds_p[t - 2];
    let c = poseidon_c();
    let s = poseidon_s();
    let m = poseidon_m();
    let p = poseidon_p();

    let component_ark = vec![vec![Fr::zero(); t]; n_rounds_f];
    let component_sigma_f = vec![vec![Fr::zero(); t]; n_rounds_f];
    let component_sigma_p = vec![Fr::zero(); n_round_p];
    let component_mix = vec![Fr::zero(); n_rounds_f - 1];
    let component_mix_s = vec![Fr::zero(); n_round_p];
    let component_mix_last = vec![Fr::zero(); n_outs];

    component_ark[0] = ark(inputs, t, c, 0);

    for j in 0..t {
        if j > 0 {
            component_ark[0][j] = inputs[j - 1];
        } else {
            component_ark[0][j] = initial_state;
        }
    }

    out
}

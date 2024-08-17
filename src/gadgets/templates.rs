// Poseidon Hash implementation from Circom Library
// https://github.com/iden3/circomlib/blob/master/circuits/poseidon.circom
// We use parameters just for ninput = 2(t = 3)! Check here for more constants
// https://github.com/iden3/circomlib/blob/master/circuits/poseidon_constants.circom

use super::*;
use halo2::halo2curves::bn256::Fr;
use utils::{poseidon_c, poseidon_m, poseidon_p, poseidon_s};

/// Exp of S-box.
pub fn sigma(item: Fr) -> Fr {
    let item_2 = item * item;
    let item_4 = item_2 * item_2;
    item * item_4
}

/// Adds round constants.
pub fn ark(state: &Vec<Fr>, t: usize, c: &Vec<Fr>, r: usize) -> Vec<Fr> {
    let mut new_state = vec![Fr::zero(); t];

    for i in 0..t {
        new_state[i] = state[i] + c[i + r]
    }

    new_state
}

/// Computes MDS matrix for MixLayer operation.
pub fn mix(state: &Vec<Fr>, t: usize, m: [[Fr; 3]; 3]) -> Vec<Fr> {
    let mut new_state = vec![Fr::zero(); t];

    for i in 0..t {
        for j in 0..t {
            new_state[i] += m[j][i] * state[j];
        }
    }

    new_state
}

/// Hash arithmetics.
pub fn poseidon_ex(inputs: [Fr; 2], initial_state: Fr, n_outs: usize) {
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

    let mut component_ark = vec![vec![Fr::zero(); t]; n_rounds_f];
    let mut component_sigma_f = vec![vec![vec![Fr::zero()]; t]; n_rounds_f];
    let component_sigma_p = vec![Fr::zero(); n_round_p];
    let mut component_mix =vec![vec![Fr::zero(); t]; n_rounds_f - 1];
    let component_mix_s = vec![Fr::zero(); n_round_p];
    let component_mix_last = vec![Fr::zero(); n_outs];

    let mut state = vec![Fr::zero(); t];

    for j in 0..t {
        if j > 0 {
            state[j] = inputs[j - 1];
        } else {
            state[j] = initial_state;
        }
    }
    
    component_ark[0] = ark(&state, t, &c, 0);        
    
    for r in 0..(n_rounds_f / 2) - 1 {        
        for j in 0..t {
            if r == 0 {
                component_sigma_f[r][j][0] = sigma(component_ark[0][j]);
            } else {
                component_sigma_f[r][j][0] = component_mix[r - 1][j];
            }           
        }               
        
        for j in 0..t {
            component_ark[r + 1][j] = component_sigma_f[r][j][0];
        }

        component_ark[r + 1] = ark(&component_ark[r + 1], t, &c, (r + 1) * t);

        for j in 0..t {
            component_mix[r][j] = component_ark[r + 1][j];
        }
        
        component_mix[r] = mix(&component_ark[r + 1], t, m);
    }   
    
    println!("{:#?}", component_mix);
}

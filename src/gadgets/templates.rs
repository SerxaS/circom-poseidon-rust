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
        new_state[i] = state[i] + c[i + r];
    }

    new_state
}

/// Computes MDS matrix for full rounds and do MixLayer operation.
pub fn mix(state: &Vec<Fr>, t: usize, m: [[Fr; 3]; 3]) -> Vec<Fr> {
    let mut new_state = vec![Fr::zero(); t];
    

    for i in 0..t {
        let mut lc = [Fr::zero()];

        for j in 0..t {
            lc[0] += m[j][i] * state[j];
        }
        
        new_state[i] = lc[0];
    }

    new_state
}

/// Computes MDS matrix for partial rounds and do MixLayer operation.
pub fn mixs(state: &Vec<Fr>, t: usize, s: &Vec<Fr>, r: usize) -> Vec<Fr> {
    let mut new_state = vec![Fr::zero(); t];

    for i in 0..t {
        new_state[0] += s[(t * 2 - 1) * r + i] * state[i];
    }

    for i in 1..t {
        new_state[i] = state[i] + state[0] * s[(t * 2 - 1) * r + t + i - 1];
    }

    new_state
}

/// Computes MDS matrix for last round and do MixLayer operation.
pub fn mixlast(state: &Vec<Fr>, t: usize, m: [[Fr; 3]; 3], s: usize) -> Fr {
    let mut new_state = vec![Fr::zero()];
    
    for j in 0..t {
        new_state[0] += m[j][s] * state[j];
    }    
    
    new_state[0]
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

    let mut component_ark = vec![vec![Fr::zero(); t]; n_rounds_f];
    let mut component_sigma_f = vec![vec![vec![Fr::zero()]; t]; n_rounds_f];
    let mut component_sigma_p = vec![Fr::zero(); n_round_p];
    let mut component_mix = vec![vec![Fr::zero(); t]; n_rounds_f - 1];
    let mut component_mix_s = vec![vec![Fr::zero(); t]; n_round_p];
    let mut component_mix_last = vec![Fr::zero(); n_outs];

    let mut state = vec![Fr::zero(); t];

    for j in 0..t {
        if j > 0 {
            state[j] = inputs[j - 1];
        } else {
            state[j] = initial_state;
        }
    }
    
    component_ark[0] = ark(&state, t, &c, 0);

    for r in 0..n_rounds_f / 2 - 1 {
        for j in 0..t {
            if r == 0 {
                component_sigma_f[r][j][0] = sigma(component_ark[0][j]);
            } else {
                component_sigma_f[r][j][0] = sigma(component_mix[r - 1][j]);
            }
        }

        let state: Vec<Fr> = component_sigma_f[r].iter().map(|item| item[0]).collect();
        component_ark[r + 1] = ark(&state, t, &c, (r + 1) * t);
        component_mix[r] = mix(&component_ark[r + 1], t, m);
    }    
    
    for j in 0..t {
        component_sigma_f[n_rounds_f / 2 - 1][j][0] = sigma(component_mix[n_rounds_f / 2 - 2][j]);
    }

    let state: Vec<Fr> = component_sigma_f[n_rounds_f / 2 - 1].iter().map(|item| item[0]).collect();
    component_ark[n_rounds_f / 2] = ark(&state, t, &c, (n_rounds_f / 2) * t);
    component_mix[n_rounds_f / 2 - 1] = mix(&component_ark[n_rounds_f / 2], t, p);
    
    for r in 0..n_round_p {
        if r == 0 {
            component_sigma_p[r] = sigma(component_mix[n_rounds_f / 2 - 1][0]);
        } else {
            component_sigma_p[r] = sigma(component_mix_s[r - 1][0]);
        }

        let mut state: Vec<Fr> = vec![Fr::zero(); t];

        for j in 0..t {
            if j == 0 {
                state[j] = component_sigma_p[r] + c[(n_rounds_f / 2 + 1) * t + r];
            } else {
                if r == 0 {
                    state[j] = component_mix[n_rounds_f / 2 - 1][j];
                } else {
                    state[j] = component_mix_s[r - 1][j];
                }
            }
        }

        component_mix_s[r] = mixs(&state, t, &s, r);
    }

    for r in 0..n_rounds_f / 2 - 1 {
        for j in 0..t {
            if r == 0 {
                component_sigma_f[n_rounds_f / 2 + r][j][0] = sigma(component_mix_s[n_round_p - 1 ][j]);
            } else {
                component_sigma_f[n_rounds_f / 2 + r][j][0] = sigma(component_mix[n_rounds_f / 2 + r - 1][j]);
            }
        }

        let state: Vec<Fr> = component_sigma_f[n_rounds_f / 2 + r].iter().map(|item| item[0]).collect();
        component_ark[n_rounds_f / 2 + r + 1] = ark(&state, t, &c, (n_rounds_f / 2 + 1) * t + n_round_p + r * t);
        component_mix[n_rounds_f / 2 + r] = mix(&component_ark[n_rounds_f / 2 + r + 1], t, m);
    }

    for j in 0..t {
        component_sigma_f[n_rounds_f - 1][j][0] = sigma(component_mix[n_rounds_f - 2][j]);
    }

    let state: Vec<Fr> = component_sigma_f[n_rounds_f - 1].iter().map(|item| item[0]).collect();
    component_mix_last[0] = mixlast(&state, t, m, 0);
        
    component_mix_last[0]
}


use std::ops::{Index, IndexMut};

use super::*;
use halo2::halo2curves::bn256::Fr;
use utils::{poseidon_c, poseidon_m, poseidon_p, poseidon_s};

/// Constructs objects.
#[derive(Clone, Debug)]
pub struct Poseidon {
    inputs: Vec<Fr>,
}

impl Poseidon {
    pub fn new(inputs: Vec<Fr>) -> Self {
        Self { inputs }
    }

    /// Exp of S-box.
    pub fn sigma(&self, item: Fr) -> Fr {
        let item_2 = item * item;
        let item_4 = item_2 * item_2;
        item * item_4
    }

    /// Adds round constants.
    pub fn ark(&self, state: &Vec<Fr>, t: usize, c: &Vec<Fr>, r: usize) -> Self {
        let mut new_state = Poseidon::new(vec![Fr::zero(); t]);

        for i in 0..t {
            new_state[i] = state[i] + c[i + r];
        }

        new_state
    }

    /// Computes MDS matrix for full rounds and do MixLayer operation.
    pub fn mix(&self, state: &Vec<Fr>, t: usize, m: [[Fr; 3]; 3]) -> Self {
        let mut new_state = Poseidon::new(vec![Fr::zero(); t]);

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
    pub fn mixs(&self, state: &Vec<Fr>, t: usize, s: &Vec<Fr>, r: usize) -> Self {
        let mut new_state = Poseidon::new(vec![Fr::zero(); t]);

        for i in 0..t {
            new_state[0] += s[(t * 2 - 1) * r + i] * state[i];
        }

        for i in 1..t {
            new_state[i] = state[i] + state[0] * s[(t * 2 - 1) * r + t + i - 1];
        }

        new_state
    }

    /// Computes MDS matrix for last round and do MixLayer operation.
    pub fn mixlast(&self, state: &Vec<Fr>, t: usize, m: [[Fr; 3]; 3], s: usize) -> Fr {
        let mut new_state = vec![Fr::zero()];

        for j in 0..t {
            new_state[0] += m[j][s] * state[j];
        }

        new_state[0]
    }

    /// Hash arithmetics.
    pub fn poseidon_ex(&self, initial_state: Fr, n_outs: usize) -> Fr {
        let n_rounds_p = [
            56, 57, 56, 60, 60, 63, 64, 63, 60, 66, 60, 65, 70, 60, 64, 68,
        ];
        let t = self.inputs.len() + 1;
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
                state[j] = self.inputs[j - 1];
            } else {
                state[j] = initial_state;
            }
        }

        component_ark[0] = self.ark(&state, t, &c, 0).inputs;

        for r in 0..n_rounds_f / 2 - 1 {
            for j in 0..t {
                if r == 0 {
                    component_sigma_f[r][j][0] = self.sigma(component_ark[0][j]);
                } else {
                    component_sigma_f[r][j][0] = self.sigma(component_mix[r - 1][j]);
                }
            }

            let state: Vec<Fr> = component_sigma_f[r].iter().map(|item| item[0]).collect();
            component_ark[r + 1] = self.ark(&state, t, &c, (r + 1) * t).inputs;
            component_mix[r] = self.mix(&component_ark[r + 1], t, m).inputs;
        }

        for j in 0..t {
            component_sigma_f[n_rounds_f / 2 - 1][j][0] =
                self.sigma(component_mix[n_rounds_f / 2 - 2][j]);
        }

        let state: Vec<Fr> = component_sigma_f[n_rounds_f / 2 - 1]
            .iter()
            .map(|item| item[0])
            .collect();
        component_ark[n_rounds_f / 2] = self.ark(&state, t, &c, (n_rounds_f / 2) * t).inputs;
        component_mix[n_rounds_f / 2 - 1] = self.mix(&component_ark[n_rounds_f / 2], t, p).inputs;

        for r in 0..n_round_p {
            if r == 0 {
                component_sigma_p[r] = self.sigma(component_mix[n_rounds_f / 2 - 1][0]);
            } else {
                component_sigma_p[r] = self.sigma(component_mix_s[r - 1][0]);
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

            component_mix_s[r] = self.mixs(&state, t, &s, r).inputs;
        }

        for r in 0..n_rounds_f / 2 - 1 {
            for j in 0..t {
                if r == 0 {
                    component_sigma_f[n_rounds_f / 2 + r][j][0] =
                        self.sigma(component_mix_s[n_round_p - 1][j]);
                } else {
                    component_sigma_f[n_rounds_f / 2 + r][j][0] =
                        self.sigma(component_mix[n_rounds_f / 2 + r - 1][j]);
                }
            }

            let state: Vec<Fr> = component_sigma_f[n_rounds_f / 2 + r]
                .iter()
                .map(|item| item[0])
                .collect();
            component_ark[n_rounds_f / 2 + r + 1] = self
                .ark(&state, t, &c, (n_rounds_f / 2 + 1) * t + n_round_p + r * t)
                .inputs;
            component_mix[n_rounds_f / 2 + r] = self
                .mix(&component_ark[n_rounds_f / 2 + r + 1], t, m)
                .inputs;
        }

        for j in 0..t {
            component_sigma_f[n_rounds_f - 1][j][0] = self.sigma(component_mix[n_rounds_f - 2][j]);
        }

        let state: Vec<Fr> = component_sigma_f[n_rounds_f - 1]
            .iter()
            .map(|item| item[0])
            .collect();
        component_mix_last[0] = self.mixlast(&state, t, m, 0);

        component_mix_last[0]
    }
}

impl Index<usize> for Poseidon {
    type Output = Fr;

    fn index(&self, idx: usize) -> &Self::Output {
        self.inputs.index(idx)
    }
}

impl IndexMut<usize> for Poseidon {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        self.inputs.index_mut(idx)
    }
}

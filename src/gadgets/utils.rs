use super::*;
use halo2::halo2curves::{bn256::Fr, ff::FromUniformBytes};

/// Returns congruent field element for the given hex string.
pub fn hex_to_field(item: &str) -> Fr {
    let item = &item[2..];
    let mut bytes = hex::decode(item).expect("Invalid parameters!");
    bytes.reverse();
    let mut temp_bytes = [0; 64];
    temp_bytes[..bytes.len()].copy_from_slice(&bytes[..]);
    Fr::from_uniform_bytes(&temp_bytes)
}

/// Returns the round constants to be used in the permutation without the prefix.
pub fn poseidon_c() -> Vec<Fr> {
    let poseidon_c_raw = constants::poseidon_c_raw();
    let poseidon_c: Vec<Fr> = poseidon_c_raw
        .iter()
        .map(|item| hex_to_field(item))
        .collect();
    poseidon_c
}

/// Returns the round constants to be used in the permutation without the prefix.
pub fn poseidon_s() -> Vec<Fr> {
    let poseidon_s_raw = constants::poseidon_s_raw();
    let poseidon_s: Vec<Fr> = poseidon_s_raw
        .iter()
        .map(|item| hex_to_field(item))
        .collect();
    poseidon_s
}

/// Returns the mds matrix elements without prefix.
pub fn poseidon_m() -> [[Fr; 3]; 3] {
    let poseidon_m = constants::poseidon_m_raw();
    poseidon_m.map(|row| row.map(|item| hex_to_field(item)))
}

/// Returns the mds matrix elements without prefix.
pub fn poseidon_p() -> [[Fr; 3]; 3] {
    let poseidon_p = constants::poseidon_p_raw();
    poseidon_p.map(|row| row.map(|item| hex_to_field(item)))
}

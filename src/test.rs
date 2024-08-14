#[cfg(test)]
mod test {
    use crate::gadgets::{templates, utils::hex_to_field};
    use halo2::halo2curves::bn256::Fr;

    #[test]
    fn test_poseidon_ark() {
        // Testing for 2 inputs.
        let inputs = [
            "0x0000000000000000000000000000000000000000000000000000000000000005",
            "0x000000000000000000000000000000000000000000000000000000000000004d",
        ]
        .map(|item| hex_to_field(item));        

        let initial_state = Fr::zero();

        templates::poseidon_ex(inputs, initial_state, 1);
    }
}

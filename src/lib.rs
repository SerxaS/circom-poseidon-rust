mod gadgets;

#[cfg(test)]
mod test {
    use crate::gadgets::{templates::Poseidon, utils::hex_to_field};
    use halo2::halo2curves::bn256::Fr;
    use num_bigint::BigInt;

    #[test]
    fn test_poseidon() {
        // Testing for 2 inputs.
        let inputs = [
            "0x0000000000000000000000000000000000000000000000000000000000000005",
            "0x000000000000000000000000000000000000000000000000000000000000004d",
        ]
        .map(|item| hex_to_field(item)).to_vec();        

        // Output taken from https://zkrepl.dev
        let expected = "6008246173323011098915936938805752727781568490715388424063708882447636047656";
        
        let poseidon = Poseidon::new(inputs);
        let initial_state = Fr::from(0);
        let n_outs = 1;        
        let output_fr = poseidon.poseidon_ex(initial_state, n_outs);
        let output_bytes = Fr::to_bytes(&output_fr);
        let output = BigInt::from_signed_bytes_le(&output_bytes).to_str_radix(10);

        assert_eq!(expected, output);        
    }
}
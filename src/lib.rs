mod gadgets;

#[cfg(test)]
mod test {
    use crate::gadgets::{templates, utils::hex_to_field};
    use halo2::halo2curves::bn256::Fr;
    use num_bigint::BigInt;

    #[test]
    fn test_poseidon() {
        // Testing for 2 inputs.
        let inputs = [
            "0x0000000000000000000000000000000000000000000000000000000000000005",
            "0x000000000000000000000000000000000000000000000000000000000000004d",
        ]
        .map(|item| hex_to_field(item));        

        // Output taken from https://zkrepl.dev
        let output = "6008246173323011098915936938805752727781568490715388424063708882447636047656";
        
        let initial_state = Fr::from(0);
        let out_fr = templates::poseidon_ex(inputs, initial_state, 1);
        let out_bytes = Fr::to_bytes(&out_fr);
        let out = BigInt::from_signed_bytes_le(&out_bytes).to_str_radix(10);

        assert_eq!(output, out);        
    }
}
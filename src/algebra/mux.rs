use algebra::prelude::*;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{prelude::*};

pub fn mux<F: Field, FG: FieldGadget<F, F>, CS: ConstraintSystem<F>>(
    mut cs: CS,
    values: &[FG],
    location: &[Boolean]) -> Result<FG, SynthesisError>
{
    let N = values.len();
    let n = location.len();
    // Assert N is a power of 2, and n = log(N)
    assert!(N & (N - 1) == 0);
    assert!(1 << n == N);
    
    let mut cur_mux_values = values.to_vec();
    for i in (0..n)
    {
        let cur_size = 1 << (n - i);
        assert!(cur_mux_values.len() == cur_size);

        let mut next_mux_values = Vec::new();
        for j in (0..cur_size).step_by(2)
        {
            let cur = FG::conditionally_select(
                cs.ns(|| format!("mux layer {:?} index {:?}", i, j)), 
                &location[n - 1 - i], 
                // true case
                &cur_mux_values[j + 1],
                // false case
                &cur_mux_values[j])?;
            next_mux_values.push(cur);
        }
        cur_mux_values = next_mux_values;
    }

    Ok(cur_mux_values[0].clone())
}

// Utility method for testing
pub fn int_to_constant_boolean_vec(
    index: u64,
    num_bits: u64) -> Vec<Boolean>
{
    let mut location = Vec::new();
    for j in (0..num_bits).rev()
    {
        location.push(Boolean::constant(index & (1 << j) != 0));
    }

    location
}

#[cfg(test)]
mod test {
    use r1cs_std::{prelude::*, test_constraint_system::TestConstraintSystem};
    use r1cs_core::{ConstraintSystem, SynthesisError};
    use crate::algebra::mux::*;

    #[test]
    fn mux_test() -> Result<(), SynthesisError> {
        use crate::alt_bn128::fr_gadget::FrGadget;
        use crate::alt_bn128::fr::Fr;
        use std::str::FromStr;

        let mut cs = TestConstraintSystem::<Fr>::new();

        let size = 8;
        let n_bits = 3;

        let mut f_gadg_vec = Vec::new();
        for i in 0..size {
            let val = Fr::from(i as u32);
            let f_gadg = FrGadget::alloc(&mut cs.ns(|| format!("generate_{:?}", i)), || Ok(val))?;
            f_gadg_vec.push(f_gadg);
        }

        for i in 0..size 
        {
            let location = int_to_constant_boolean_vec(i, n_bits);
            let mux_res = mux(cs.ns(|| format!("mux {:?}", i)), &f_gadg_vec, &location)?;
            mux_res.enforce_equal(&mut cs.ns(|| format!("test {:?}", i)), &f_gadg_vec[i as usize]);
        }
        if !cs.is_satisfied() {
            println!("{:?}", cs.which_is_unsatisfied().unwrap());
        }

        println!("{:?}", cs.num_constraints());
        cs.print_named_objects();
        assert!(cs.is_satisfied());
        
        println!("{:?}", cs.num_constraints());

        Ok(())
    }
}

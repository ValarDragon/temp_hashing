use algebra::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;

/// Struct describing polynomials
/// Consider also having it take a sparse polynomial of fixed coefficients
#[derive(Clone)]
pub struct DensePolynomial<F: Field, FG: FieldGadget<F, F>>
{
    pub coeffs: Vec<FG>,
    pub _phantom: F,
}

impl<F : Field, FG: FieldGadget<F, F>> DensePolynomial<F, FG> {
    /// Evaluates the constraints and just gives you the gadget for the result.
    /// Caution for use in holographic lincheck: The output has 2 entries in one matrix
    pub fn evaluate_constraints<CS: ConstraintSystem<F>>(&self, 
        mut cs: CS,
        x: &FG) -> Result<FG, SynthesisError>
        where FG: FieldGadget<F, F>
    {
        let mut p_cs = &mut cs.ns(|| "polynomial eval");
        let mut res = FG::zero(&mut p_cs)?;
        let mut cur_pow_x = FG::one(&mut p_cs)?;
        for i in 0..self.coeffs.len()
        {
            let term = cur_pow_x.mul(p_cs.ns(|| format!("c[i] * x^{:?}", i)), &self.coeffs[i])?;
            res.add_in_place(&mut p_cs, &term);

            cur_pow_x = cur_pow_x.mul(p_cs.ns(|| format!("x^{:?}", i)), x)?;
        }
        Ok(res)
    }
}

#[test]
fn test_polynomial_constraints() {
    use r1cs_std::alloc::*;
    
    use r1cs_std::eq::EqGadget;
    use r1cs_std::{test_constraint_system::TestConstraintSystem};
    use crate::alt_bn128::fr_gadget::FrGadget;
    use crate::alt_bn128::fr::Fr;
    use std::str::FromStr;
    use algebra::prelude::*;
    
    use crate::algebra::polynomial::DensePolynomial;

    let mut cs = TestConstraintSystem::<Fr>::new();

    // 1 + 0x + 1 x^2
    let coeffs = vec![FrGadget::one(&mut cs).unwrap(), FrGadget::zero(&mut cs).unwrap(), FrGadget::one(&mut cs).unwrap()];
    let poly = DensePolynomial{
        coeffs,
        _phantom: Fr::zero(),
    };
    let x_fr: Fr = Fr::from(5u32);
    let exp_fr = Fr::from(26u32);

    let x = FrGadget::alloc(&mut cs.ns(|| "generate_x"), || Ok(x_fr)).unwrap();
    let exp = FrGadget::alloc(&mut cs.ns(|| "generate_exp"), || Ok(exp_fr)).unwrap();
    
    let res = poly.evaluate_constraints(&mut cs, &x).unwrap();
    
    res.enforce_equal(&mut cs, &exp).unwrap();
    if !cs.is_satisfied() {
        println!("{:?}", cs.which_is_unsatisfied().unwrap());
    }
    println!("{:?}", cs.num_constraints());
    cs.print_named_objects();
    assert!(cs.is_satisfied());
}

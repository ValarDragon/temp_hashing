use algebra::Field;
use r1cs_core::ConstraintSystem;
use r1cs_std::prelude::*;

/// Struct describing vanishing polynomials for a multiplicative coset H,
/// with |H| a power of 2.
/// As H is a coset, every element can be described as h*g^i,
/// and therefore has vanishing polynomial Z_H(x) = x^|H| - h^|H|
pub struct VanishingPolynomial<F: Field>
{
    /// h^|H|
    constant_term: F,
    /// log_2(|H|)
    dim_h: u64,
    // |H|
    order_h: u64,
}

impl<F : Field> VanishingPolynomial<F> {
    pub fn new(offset: F, dim_h: u64) -> Self
    {
        let order_h = 1 << dim_h;
        let vp = VanishingPolynomial{
            constant_term: offset.pow([order_h]),
            dim_h,
            order_h,
        };
        vp
    }

    pub fn evaluate(&self, x: &F) -> F {
        let mut result = x.pow([self.order_h]);
        result -= &self.constant_term;
        result
    }

    /// Evaluates the constraints and just gives you the gadget for the result.
    /// Caution for use in holographic lincheck: The output has 2 entries in one matrix
    pub fn evaluate_constraints<CS: ConstraintSystem<F>, FG>(&self, 
        mut cs: CS,
        x: &FG) -> FG
        where FG: FieldGadget<F, F>
    {
        let mut vp_cs = &mut cs.ns(|| "vanishing polynomial");
        if self.dim_h == 1
        {
            let result = x.sub(&mut vp_cs.ns(|| "compute result"), x).unwrap();
            return result;
        }
        let mut cur = x.square(vp_cs.ns(|| format!("compute x^(2^{:?})", 1))).unwrap();
        for i in 1..self.dim_h 
        {
            cur.square_in_place(vp_cs.ns(|| format!("compute x^(2^{:?})", i + 1))).unwrap();
        }
        cur.sub_constant_in_place(vp_cs.ns(|| "compute result"), &self.constant_term).unwrap();
        cur
    }
}

#[test]
fn test_vanishing_polynomial() {
    use crate::alt_bn128::fr::Fr;
    use crate::algebra::vanishing_polynomial::VanishingPolynomial;
    use std::str::FromStr;
    

    let dim_h = 30;
    let h: Fr = Fr::from(372u32);
    let vp = VanishingPolynomial::new(h, dim_h);
    let x: Fr = Fr::from(17u32);
    let res = vp.evaluate(&x);
    
    // evaluated in sage
    let exp = Fr::from_str("13050942103627446995176210260937927095935575529594931712231973977461351072286").map_err(|_| ()).unwrap();

    assert_eq!(res, exp);
}

#[test]
fn test_vanishing_polynomial_constraints() {
    use r1cs_std::alloc::*;
    
    use r1cs_std::eq::EqGadget;
    use r1cs_std::{test_constraint_system::TestConstraintSystem};
    use crate::alt_bn128::fr_gadget::FrGadget;
    use crate::alt_bn128::fr::Fr;
    use std::str::FromStr;
    
    use crate::algebra::vanishing_polynomial::VanishingPolynomial;

    let mut cs = TestConstraintSystem::<Fr>::new();

    let dim_h = 30;
    let h: Fr = Fr::from(372u32);
    let vp = VanishingPolynomial::new(h, dim_h);
    let x_fr: Fr = Fr::from(17u32);
    let x = FrGadget::alloc(&mut cs.ns(|| "generate_x"), || Ok(x_fr)).unwrap();
    let res = vp.evaluate_constraints(&mut cs, &x);
    
    // evaluated in sage
    let exp_fr = Fr::from_str("13050942103627446995176210260937927095935575529594931712231973977461351072286").map_err(|_| ()).unwrap();
    let exp = FrGadget::alloc(&mut cs.ns(|| "generate_exp"), || Ok(exp_fr)).unwrap();

    res.enforce_equal(&mut cs, &exp).unwrap();
    if !cs.is_satisfied() {
        println!("{:?}", cs.which_is_unsatisfied().unwrap());
    }
    println!("{:?}", cs.num_constraints());
    cs.print_named_objects();
    assert!(cs.is_satisfied());
    // library does not optimize away the final enforce_equal constraint
    assert!(cs.num_constraints() == (dim_h as usize) + 1);
}

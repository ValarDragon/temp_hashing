use algebra::{Field, PrimeField};
use r1cs_core::ConstraintSystem;
use r1cs_std::prelude::*;
use r1cs_std::fields::fp::FpGadget;
use crate::algebra::vanishing_polynomial::VanishingPolynomial;
use algebra::fields::batch_inversion;
use num_traits::Zero;
use std::convert::{AsRef, From};

/// Struct describing Lagrange interpolation for a multiplicative coset I,
/// with |I| a power of 2.
/// TODO: Pull in lagrange poly explanation from libiop
pub struct LagrangeInterpolator<F: PrimeField>
{
    domain_order : usize,
    all_domain_elems: Vec<F>,
    v_inv_elems: Vec<F>,
    domain_vp: VanishingPolynomial<F>,
    poly_evaluations: Vec<F>,
}

/// TODO: My notes on how to do the constraints.
/// We can construct the lagrange coefficients in just k constraints!
/// For coeff[i], its `(v_inv[i] * x - v_inv[i] * elem[i]) * (coeff[i]) = Z_I(x)
/// (We compute Z_I(x) re-using coefficients of Z_H(x), so it only takes 1 more addition when optimized.)
/// Then we use k constraints for summing & multiplying these with the values of f.
/// 
/// Struct describing a Lagrange interpolation gadget for a multiplicative coset I,
/// with |I| a power of 2.
/// TODO: Pull in lagrange poly explanation from libiop
pub struct LagrangeInterpolationGadget<F: PrimeField>
{
    pub lagrange_interpolator : LagrangeInterpolator<F>,
    // A hack for optimizations outside of the lagrange interpolation
    pub vp_t : Option<FpGadget<F>>,
    poly_evaluations: Vec<FpGadget<F>>,
}

impl<F: PrimeField> LagrangeInterpolator<F> {
    pub fn new(
        domain_offset: F,
        domain_generator: F, 
        domain_dim: u64, 
        poly_evaluations: Vec<F>) -> Self
    {
        let domain_order = 1 << domain_dim;
        assert_eq!(poly_evaluations.len(), domain_order);
        let mut cur_elem = domain_offset;
        let mut all_domain_elems = vec![domain_offset];
        let mut v_inv_elems : Vec<F> = Vec::new();
        /// Cache all elements in the domain
        for _ in 1..domain_order
        {
            cur_elem *= domain_generator;
            all_domain_elems.push(cur_elem);
        }
        /// By computing the following elements as constants, 
        /// we can further reduce the interpolation costs.
        /// 
        /// m = order of the interpolation domain
        /// v_inv[i] = prod_{j != i} h(g^i - g^j)
        /// We use the following facts to compute this:
        ///   v_inv[0] = m*h^{m-1}
        ///   v_inv[i] = g^{-1} * v_inv[i-1]
        /// TODO: Include proof of the above two points
        let g_inv = domain_generator.inverse().unwrap();
        let m = F::from((1 << domain_dim) as u64);
        let mut v_inv_i = m * domain_offset.pow([(domain_order - 1) as u64]);
        for _ in 0..domain_order
        {
            v_inv_elems.push(v_inv_i);
            v_inv_i *= g_inv;
        }

        /// TODO: Ideally we'd cache the intermediate terms with Z_H(x) evaluations, since most of the exponents are precomputed.
        let vp = VanishingPolynomial::new(domain_offset, domain_dim);

        let lagrange_interpolation : LagrangeInterpolator<F> = LagrangeInterpolator{
            domain_order,
            all_domain_elems,
            v_inv_elems,
            domain_vp: vp,
            poly_evaluations,
        };
        lagrange_interpolation
    }

    fn compute_lagrange_coefficients(&self, interpolation_point: F) -> Vec<F>
    {
        /*
        * Let t be the interpolation point, H be the multiplicative coset, with elements of the form h*g^i.
        Compute each L_{i,H}(t) as Z_{H}(t) * v_i / (t- h g^i)
        where:
        - Z_{H}(t) = \prod_{j} (t-h*g^j) = (t^m-h^m), and
        - v_{i} = 1 / \prod_{j \neq i} h(g^i-g^j).
        Below we use the fact that v_{0} = 1/(m * h^(m-1)) and v_{i+1} = g * v_{i}.
        We compute the inverse of each coefficient, and then batch invert the entire result.
        TODO: explain deriviation more step by step
        */
        // TODO: Implement batch_inverse & mul like libiop for better efficiency
        let vp_t_inv = self.domain_vp.evaluate(&interpolation_point).inverse().unwrap();
        let mut inverted_lagrange_coeffs : Vec<F> = Vec::with_capacity(self.all_domain_elems.len());
        for i in 0..self.domain_order
        {
            let l = vp_t_inv * self.v_inv_elems[i];
            let r = self.all_domain_elems[i];
            inverted_lagrange_coeffs.push(l * (interpolation_point - r));
        }
        let lagrange_coeffs = inverted_lagrange_coeffs.as_mut_slice();
        batch_inversion::<F>(lagrange_coeffs);
        lagrange_coeffs.iter().cloned().collect()
    }

    pub fn interpolate(&self, interpolation_point: F) -> F
    {
        let lagrange_coeffs = self.compute_lagrange_coefficients(interpolation_point);
        let mut interpolation = F::zero();
        for i in 0..self.domain_order
        {
            interpolation += (lagrange_coeffs[i] * self.poly_evaluations[i]);
        }
        interpolation
    }
}

impl<F: PrimeField> LagrangeInterpolationGadget<F> {
    pub fn new(
        domain_offset: F, 
        domain_generator: F, 
        domain_dim: u64, 
        poly_evaluations: Vec<FpGadget<F>>) -> Self
    {
        let mut poly_evaluations_F : Vec<F> = Vec::new();
        for i in 0..(1 << domain_dim)
        {
            poly_evaluations_F.push(poly_evaluations[i].get_value().unwrap());
        }

        let lagrange_interpolator : LagrangeInterpolator<F> = LagrangeInterpolator::new(
            domain_offset,
            domain_generator,
            domain_dim,
            poly_evaluations_F,
        );

        let lagrange_interpolation_gadget = LagrangeInterpolationGadget{
            lagrange_interpolator,
            vp_t : None,
            poly_evaluations,
        };
        lagrange_interpolation_gadget
    }

    fn compute_lagrange_coefficients_constraints<CS: ConstraintSystem<F>>(&mut self,
        mut cs: CS,
        interpolation_point: &FpGadget<F>,
        ) -> Vec<FpGadget<F>>
    {
        let t = interpolation_point;
        let lagrange_coeffs = self.lagrange_interpolator.compute_lagrange_coefficients(t.get_value().unwrap());
        let mut lagrange_coeffs_FG : Vec<FpGadget<F>> = Vec::new();
        // Now we convert these lagrange coefficients to gadgets, and then constrain them.
        // The i-th lagrange coefficients constraint is:
        // (v_inv[i] * t - v_inv[i] * domain_elem[i]) * (coeff) = 1/Z_I(t)
        // 
        let vp_t = self.lagrange_interpolator.domain_vp.evaluate_constraints(&mut cs, t);
        let inv_vp_t = vp_t.inverse(cs.ns(|| "Take inverse of Z_I(t)")).unwrap();
        self.vp_t = Some(vp_t);
        for i in 0..(self.lagrange_interpolator.domain_order)
        {
            let constant = (-self.lagrange_interpolator.all_domain_elems[i]) * self.lagrange_interpolator.v_inv_elems[i];
            let mut A_element = t.mul_by_constant(&mut cs, &self.lagrange_interpolator.v_inv_elems[i]).unwrap();
            A_element.add_constant_in_place(&mut cs, &constant);
            
            let lag_coeff = FpGadget::<F>::alloc(
                &mut cs.ns(|| format!("generate lagrange coefficient {:?}", i)), 
                || Ok(lagrange_coeffs[i])).unwrap();
                lagrange_coeffs_FG.push(lag_coeff);
            // Enforce the actual constraint (A_element) * (lagrange_coeff) = 1/Z_I(t)
            A_element.mul_equals(cs.ns(|| format!("Check the {:?}th lagrange coefficient", i)), 
                &lagrange_coeffs_FG[i], 
                &inv_vp_t).unwrap();
        }
        return lagrange_coeffs_FG;
    }

    pub fn interpolate_constraints<CS: ConstraintSystem<F>>(&mut self, 
        mut cs: CS, 
        interpolation_point: &FpGadget<F>) -> FpGadget<F>
    {
        let lagrange_coeffs = self.compute_lagrange_coefficients_constraints(&mut cs, interpolation_point);
        let mut interpolation = FpGadget::<F>::from(&mut cs, &F::zero());
        // Set interpolation to be: sum_{i in domain} lagrange_coeff(i)*f(i)
        for i in 0..self.lagrange_interpolator.domain_order
        {
            let intermediate = lagrange_coeffs[i].mul(
                cs.ns(|| format!("Compute the product of {:?}th lagrange coefficient and polynomial interpolation", i)),
                &self.poly_evaluations[i]).unwrap();
            interpolation = interpolation.add(&mut cs, &intermediate).unwrap();
        }
        interpolation
    }
}

#[cfg(test)]
mod test {
    use r1cs_std::{prelude::*, test_constraint_system::TestConstraintSystem};
    use r1cs_core::ConstraintSystem;
    use crate::algebra::lagrange_interpolation::*;
    use crate::algebra::domain::Domain;
    use algebra_core::fields::*;

    use r1cs_std::alloc::*;
    use r1cs_std::eq::EqGadget;
    use crate::alt_bn128::fr_gadget::FrGadget;
    use crate::alt_bn128::fr::Fr;
    use std::str::FromStr;
    use num_traits::*;

    #[test]
    fn single_round_fri_test() {
        use crate::algebra::polynomial::DensePolynomial;

        let mut cs = TestConstraintSystem::<Fr>::new();

        // 1 + 0x + 0x^2 + 1 x^3
        let coeffs = vec![FrGadget::one(&mut cs).unwrap(), 
            FrGadget::zero(&mut cs).unwrap(), 
            FrGadget::zero(&mut cs).unwrap(), 
            FrGadget::one(&mut cs).unwrap()];
        let poly = DensePolynomial::<Fr, FrGadget>{   
            coeffs,
            _phantom: Fr::zero(),
        };

        let root = Fr::two_adic_root_of_unity();
        let subgroup_generator = root.pow(&[1u64 << 26]);
        assert_eq!(root.pow(&[1 << 28]), Fr::one());
        assert_eq!(subgroup_generator.pow(&[4]), Fr::one());
        let subgroup = vec![Fr::one(), subgroup_generator, subgroup_generator.pow(&[2]), subgroup_generator.pow(&[3])];

        let mut oracle_evals = Vec::new();
        let domain_offset = FrGadget::alloc(cs.ns(|| "domain offset"), || Ok(Fr::multiplicative_generator())).unwrap();
        for i in 0..4
        {
            let mut cs_i = cs.ns(|| format!("{:?}", i));
            let coset_point = domain_offset.mul_by_constant(&mut cs_i, &subgroup_generator.pow(&[i])).unwrap();
            oracle_evals.push(poly.evaluate_constraints(&mut cs_i, &coset_point).unwrap());
        }
        let alpha = Fr::from(11u64);
        let alpha_FG = FrGadget::alloc(&mut cs.ns(||"alloc alpha"), || Ok(alpha)).unwrap();
        let expected = poly.evaluate_constraints(&mut cs, &alpha_FG).unwrap().get_value().unwrap();

        let L = Domain{
            gen: subgroup_generator,
            offset: Fr::multiplicative_generator(),
            dim: 2,
        };
        // Check the single round FRI prover as well
        let interpolator = LagrangeInterpolator::new(
            L.offset,
            L.gen,
            2,
            oracle_evals.into_iter().map(|x| x.get_value().unwrap()).collect());
        assert_eq!(expected, interpolator.interpolate(alpha));
    }
}
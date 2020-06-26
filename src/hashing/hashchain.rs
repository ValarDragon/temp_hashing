use algebra::prelude::*;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{prelude::*, fields::fp::FpGadget};

use crate::hashing::*;
use crate::hashing::sponge::*;
use num_traits::*;

pub struct HashChain<F: PrimeField, P: Permutation<F> + Clone>
{
    sponge: AlgebraicSponge<F, P>,
}

impl<F: PrimeField, P: Permutation<F> + Clone> HashChain<F, P>
{
    pub fn new(sponge: AlgebraicSponge<F, P>) -> Self
    {
        let mut sponge_c = sponge.clone();
        sponge_c.reset();
        HashChain{sponge: sponge_c}
    }

    pub fn absorb(&mut self, elems: &[F])
    {
        self.sponge.absorb(elems);
    }

    pub fn squeeze(&mut self, num_elements: usize) -> Vec<F>
    {
        self.sponge.squeeze(num_elements)
    }

    pub fn squeeze_ints(&mut self, num_elements: usize, bits_per_elem: usize) -> Vec<u64>
    {
        let squeezed_field_elems = self.sponge.squeeze(num_elements);

        let mut vec_of_bits : Vec<Vec<bool>> = squeezed_field_elems.into_iter().map(|x| x.into_repr().to_bits()).collect();
        print!("{:?}", vec_of_bits);
        let mut squeezed_ints = Vec::new();
        for i in 0..num_elements {
            let mut cur = 0;
            vec_of_bits[i].reverse();
            for j in 0..bits_per_elem {
                cur <<= 1;
                cur += vec_of_bits[i][j] as u64;
            }
            squeezed_ints.push(cur);
        }

        squeezed_ints
    }
}

pub struct HashChainGadget<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG> + Clone>
{
    sponge: AlgebraicSpongeGadget<F, FG, P>,
}

impl<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG> + Clone> HashChainGadget<F, FG, P>
{
    pub fn new<CS: ConstraintSystem<F>>(mut cs: CS, sponge: AlgebraicSpongeGadget<F, FG, P>) -> Self
    {
        let mut sponge_c = sponge.clone();
        sponge_c.reset(&mut cs);
        HashChainGadget{sponge: sponge_c}
    }

    pub fn absorb<CS: ConstraintSystem<F>>(&mut self, mut cs: CS, elems: &[FG])
    {
        self.sponge.absorb(&mut cs, elems);
    }

    pub fn squeeze<CS: ConstraintSystem<F>>(&mut self, mut cs: CS, num_elements: usize) -> Result<Vec<FG>, SynthesisError>
    {
        self.sponge.squeeze(&mut cs, num_elements)
    }

    pub fn squeeze_ints<CS: ConstraintSystem<F>>(&mut self, mut cs: CS, num_elements: usize, bits_per_elem: usize) -> Result<Vec<Vec<Boolean>>, SynthesisError>
    {
        let mut cs_cur = cs.ns(|| "Squeeze ints");
        let squeezed_field_elems = self.sponge.squeeze(&mut cs_cur, num_elements)?;

        let mut vec_of_bits = Vec::<Vec<Boolean>>::new();
        for i in 0..squeezed_field_elems.len(){
            vec_of_bits.push(squeezed_field_elems[i].to_bits(&mut cs_cur.ns(|| format!("field elem {:?}", i))).unwrap());
        }
        let elem_len = vec_of_bits[0].len();
        for i in 0..num_elements {
            let mut cur = 0;
            vec_of_bits[i] = vec_of_bits[i][0..bits_per_elem].to_vec();
        }

        Ok(vec_of_bits)
    }
}


#[cfg(test)]
mod test {
    use r1cs_std::{prelude::*, test_constraint_system::TestConstraintSystem};
    use r1cs_core::ConstraintSystem;
    use crate::algebra::{ mux::*};

    use algebra::{FpParameters, prelude::*};
    use r1cs_std::prelude::*;
    use r1cs_core::SynthesisError;
    use r1cs_std::alloc::*;
    use r1cs_std::eq::EqGadget;
    use crate::alt_bn128::fr_gadget::FrGadget;
    use crate::alt_bn128::fr::Fr;
    use std::str::FromStr;
    use crate::hashing::{*, dummy_permutation::*, sponge::*, hashchain::*};

    #[test]
    fn hashchain_consistency_test() -> Result<(), SynthesisError> {
        let rate = 3;
        let capacity = 1;
        let P = DummyPermutation{};
        let sponge = AlgebraicSponge::<Fr, DummyPermutation>::new(rate, capacity, P);
        let mut hashchain = HashChain::new(sponge);
        hashchain.absorb(&[Fr::one(), Fr::one() + Fr::one()]);
        let res = hashchain.squeeze_ints(2, 1);
        assert_eq!(res, vec![0, 1]);        
        Ok(())
    }
}
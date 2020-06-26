use algebra::prelude::*;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{prelude::*, fields::fp::FpGadget};

use crate::hashing::*;
use crate::hashing::sponge::*;
use num_traits::*;

pub struct TwoToOneHash<F: PrimeField, P: Permutation<F> + Clone>
{
    sponge: AlgebraicSponge<F, P>,
}

impl<F: PrimeField, P: Permutation<F> + Clone> TwoToOneHash<F, P>
{
    pub fn new(sponge: AlgebraicSponge<F, P>) -> Self
    {
        let mut sponge_c = sponge.clone();
        sponge_c.reset();
        TwoToOneHash{sponge: sponge_c}
    }

    pub fn hash(&self, left: F, right : F) -> F
    {
        let mut sponge_c = self.sponge.clone();
        sponge_c.absorb(&[left, right]);
        // TODO: Make generic for smaller fields
        sponge_c.squeeze(1)[0]
    }
}

pub struct TwoToOneHashGadget<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG> + Clone>
{
    sponge: AlgebraicSpongeGadget<F, FG, P>,
}

impl<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG> + Clone> TwoToOneHashGadget<F, FG, P>
{
    pub fn new<CS: ConstraintSystem<F>>(mut cs: CS, sponge: AlgebraicSpongeGadget<F, FG, P>) -> Self
    {
        let mut sponge_c = sponge.clone();
        sponge_c.reset(&mut cs);
        TwoToOneHashGadget{sponge: sponge_c}
    }

    pub fn hash<CS: ConstraintSystem<F>>(&self, mut cs: CS, left: FG, right : FG) -> Result<FG, SynthesisError>
    {
        let mut sponge_c = self.sponge.clone();
        sponge_c.absorb(&mut cs, &[left, right]);
        // TODO: Make generic for smaller fields
        Ok(sponge_c.squeeze(&mut cs, 1)?[0].clone())
    }
}
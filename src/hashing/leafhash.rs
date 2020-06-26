use algebra::prelude::*;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{prelude::*, fields::fp::FpGadget};

use crate::hashing::*;
use crate::hashing::sponge::*;
use num_traits::*;

pub struct LeafHash<F: PrimeField, P: Permutation<F> + Clone>
{
    sponge: AlgebraicSponge<F, P>,
}

impl<F: PrimeField, P: Permutation<F> + Clone> LeafHash<F, P>
{
    pub fn new(sponge: AlgebraicSponge<F, P>) -> Self
    {
        let mut sponge_copy = sponge.clone();
        sponge_copy.reset();
        LeafHash{sponge: sponge_copy}
    }

    pub fn hash(&self, elems: &[F]) -> F
    {
        let mut sponge_copy = self.sponge.clone();
        sponge_copy.absorb(elems);
        // TODO: Make generic for smaller fields
        sponge_copy.squeeze(1)[0]
    }

    pub fn zk_hash(&self, elems: &[F], salt: F) -> F
    {
        let mut sponge_copy = self.sponge.clone();
        sponge_copy.absorb(elems);
        sponge_copy.absorb(&[salt]);
        // TODO: Make generic for smaller fields
        sponge_copy.squeeze(1)[0]
    }
}

pub struct LeafHashGadget<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG> + Clone>
{
    sponge: AlgebraicSpongeGadget<F, FG, P>,
}

impl<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG> + Clone> LeafHashGadget<F, FG, P>
{
    pub fn new<CS: ConstraintSystem<F>>(mut cs: CS, sponge: AlgebraicSpongeGadget<F, FG, P>) -> Self
    {
        let mut sponge_copy = sponge.clone();
        sponge_copy.reset(&mut cs);
        LeafHashGadget{sponge: sponge_copy}
    }

    pub fn hash<CS: ConstraintSystem<F>>(&self, mut cs: CS, elems: &[FG]) -> Result<FG, SynthesisError>
    {
        let mut sponge_copy = self.sponge.clone();
        sponge_copy.absorb(&mut cs, elems);
        // TODO: Make generic for smaller fields
        Ok(sponge_copy.squeeze(&mut cs, 1)?[0].clone())
    }

    pub fn zk_hash<CS: ConstraintSystem<F>>(&self, mut cs: CS, elems: &[FG], salt: FG) -> Result<FG, SynthesisError>
    {
        let mut sponge_copy = self.sponge.clone();
        sponge_copy.absorb(&mut cs, elems);
        sponge_copy.absorb(&mut cs, &[salt]);
        // TODO: Make generic for smaller fields
        Ok(sponge_copy.squeeze(&mut cs, 1)?[0].clone())
    }
}
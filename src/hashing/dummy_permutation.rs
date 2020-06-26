use algebra::fields::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{prelude::*};
use crate::hashing::*;

#[derive(Clone)]
pub struct DummyPermutation {}

impl<F: Field> Permutation<F> for DummyPermutation
{
    fn permute(&self, state: &mut[F])
    {
        for i in 0..state.len()
        {
            state[i] += F::one();
        }
    }
}

#[derive(Clone)]
pub struct SeededDummyPermutation<F: Field> {pub seed: F}

impl<F: Field> Permutation<F> for SeededDummyPermutation<F>
{
    fn permute(&self, state: &mut[F])
    {
        let mut cur = self.seed;
        for i in 0..state.len()
        {
            state[i] += cur;
            cur += self.seed;
        }
    }
}

impl<F: Field, FG: FieldGadget<F,F>> PermutationGadget<F, FG> for SeededDummyPermutation<F>
{
    fn permute<CS: ConstraintSystem<F>>(&self, mut cs: CS, state: &mut[FG]) -> Result<(), SynthesisError>
    {
        let mut cur = self.seed;
        for i in 0..state.len()
        {
            state[i].add_constant_in_place(&mut cs, &cur)?;
            cur += self.seed;
        }
        Ok(())
    }
}

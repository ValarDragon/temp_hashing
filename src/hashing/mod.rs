pub mod sponge;
pub mod leafhash;
pub mod hashchain;
pub mod poseidon;
pub mod rescue;
pub mod dummy_permutation;
pub mod two_to_one_hash;

use algebra::fields::Field;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{prelude::*};

// Const Generics aren't stable, so the size parameter cannot be templated.
pub trait Permutation<F: Field>
{
    fn permute(
        &self,
        state: &mut [F]);

    // pub fn print_soundness(&self);
}

pub trait PermutationGadget<F, FG>
where
    F: Field, FG: FieldGadget<F, F>
{
    fn permute<CS: ConstraintSystem<F>>(
        &self, 
        cs: CS, 
        state: &mut [FG]) -> Result<(), SynthesisError>;
}
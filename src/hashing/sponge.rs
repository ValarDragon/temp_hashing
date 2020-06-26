use algebra::prelude::*;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{prelude::*, fields::fp::FpGadget};

use crate::hashing::*;
use num_traits::*;

#[derive(Clone)]
enum SpongeState {
    Absorbing{next_absorb_index : usize},
    Squeezing{next_squeeze_index : usize},
}

#[derive(Clone)]
pub struct AlgebraicSponge<F: PrimeField, P: Permutation<F>>
{
    state: Vec<F>,
    rate: usize,
    capacity: usize,
    permutation: P,
    mode: SpongeState,
}

impl<F: PrimeField, P: Permutation<F>> AlgebraicSponge<F, P>
{
    pub fn new(rate: usize, capacity: usize, permutation: P) -> Self
    {
        let state = vec![F::zero(); rate + capacity];
        let mode = SpongeState::Absorbing{ next_absorb_index: 0};

        AlgebraicSponge{
            state,
            rate,
            capacity,
            permutation,
            mode,
        }
    }

    pub fn reset(&mut self)
    {
        self.state = vec![F::zero(); self.rate + self.capacity];
        self.mode = SpongeState::Absorbing{ next_absorb_index: 0};
    }

    fn permute(&mut self)
    {
        self.permutation.permute(&mut self.state);
    }

    pub fn absorb(&mut self, elements: &[F])
    {
        match self.mode {
            SpongeState::Absorbing{ next_absorb_index } => {
                let mut absorb_index = next_absorb_index;
                if (absorb_index == self.rate)
                {
                    self.permute();
                    absorb_index = 0;
                }
                self.absorb_internal(absorb_index, elements);
            },
            SpongeState::Squeezing{ next_squeeze_index } => { 
                self.permute(); 
                self.absorb_internal(0, elements);
            },
        };
    }

    // Absorbs everything in elements, this does not end in an absorbtion.
    fn absorb_internal(&mut self, rate_start_index: usize, elements: &[F])
    {
        // if we can finish in this call
        if rate_start_index + elements.len() <= self.rate 
        {
            for i in 0..elements.len()
            {
                self.state[i + rate_start_index] += elements[i];
            }
            self.mode = SpongeState::Absorbing{ next_absorb_index: rate_start_index + elements.len() };
            return;
        }
        // otherwise absorb (rate - rate_start_index) elements
        let num_elements_absorbed = self.rate - rate_start_index;
        for i in 0..num_elements_absorbed
        {
            self.state[i + rate_start_index] += elements[i];
        }
        self.permute();
        // Tail recurse, with the input elements being truncated by num elements absorbed
        self.absorb_internal(0, &elements[num_elements_absorbed..]);
    }

    pub fn squeeze(&mut self, num_elements: usize) -> Vec<F>
    {
        let mut squeezed_elems = vec![F::zero(); num_elements];
        match self.mode {
            SpongeState::Absorbing{ next_absorb_index } => { 
                self.permute();
                self.squeeze_internal(0, &mut squeezed_elems);
            },
            SpongeState::Squeezing{ next_squeeze_index } => {
                let mut squeeze_index = next_squeeze_index;
                if (squeeze_index == self.rate)
                {
                    self.permute();
                    squeeze_index = 0;
                }
                self.squeeze_internal(squeeze_index, &mut squeezed_elems);
            },
        };

        squeezed_elems
    }

    // Squeeze |output| many elements. This does not end in a squeeze
    fn squeeze_internal(&mut self, rate_start_index: usize, output: &mut[F])
    {
        // if we can finish in this call
        if rate_start_index + output.len() <= self.rate 
        {
            for i in 0..output.len()
            {
                output[i] = self.state[i + rate_start_index];
            }
            self.mode = SpongeState::Squeezing{ next_squeeze_index: rate_start_index + output.len() };
            return;
        }
        // otherwise squeeze (rate - rate_start_index) elements
        let num_elements_squeezed = self.rate - rate_start_index;
        for i in 0..num_elements_squeezed
        {
            output[i] = self.state[i + rate_start_index];
        }

        // Unless we are done with squeezing in this call, permute.
        if (output.len() != self.rate)
        {
            self.permute();
        }
        // Tail recurse, with the correct change to indices in output happening due to changing the slice
        self.squeeze_internal(0, &mut output[num_elements_squeezed..]);
    }
}

#[derive(Clone)]
pub struct AlgebraicSpongeGadget<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG>>
{
    state: Vec<FG>,
    rate: usize,
    capacity: usize,
    permutation: P,
    mode: SpongeState,

    // needed for cs unique prefixes
    global_permute_index: usize,
    _phantom: F,
}

impl<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG>> AlgebraicSpongeGadget<F, FG, P>
{
    fn zero_state<CS: ConstraintSystem<F>>(mut cs: CS, state_size: usize) -> Result<Vec<FG>, SynthesisError>
    {
        // Initialize state to a vector of 0
        let mut state = Vec::new();
        for i in 0..state_size
        {
            state.push(FG::zero(&mut cs)?);
        }

        Ok(state)
    }
    pub fn new<CS: ConstraintSystem<F>>(mut cs: CS, rate: usize, capacity: usize, permutation: P) -> Result<Self, SynthesisError>
    {
        let state = AlgebraicSpongeGadget::<F, FG, P>::zero_state(&mut cs, rate + capacity)?;
        let mode = SpongeState::Absorbing{ next_absorb_index: 0};

        Ok(AlgebraicSpongeGadget{
            state,
            rate,
            capacity,
            permutation,
            mode,
            global_permute_index: 0,
            _phantom: F::zero(),
        })
    }

    pub fn reset<CS: ConstraintSystem<F>>(&mut self, mut cs: CS) -> Result<(), SynthesisError>
    {
        self.state = AlgebraicSpongeGadget::<F, FG, P>::zero_state(&mut cs, self.rate + self.capacity)?;
        self.mode = SpongeState::Absorbing{ next_absorb_index: 0};
        Ok(())
    }

    fn permute<CS: ConstraintSystem<F>>(&mut self, mut cs: CS) -> Result<(), SynthesisError>
    {
        self.permutation.permute(&mut cs.ns(|| format!("permute {:?}", self.global_permute_index)), &mut self.state)?;
        self.global_permute_index += 1;
        Ok(())
    }

    pub fn absorb<CS: ConstraintSystem<F>>(&mut self, mut cs: CS, elements: &[FG]) -> Result<(), SynthesisError>
    {
        match self.mode {
            SpongeState::Absorbing{ next_absorb_index } => {
                let mut absorb_index = next_absorb_index;
                if (absorb_index == self.rate)
                {
                    self.permute(&mut cs)?;
                    absorb_index = 0;
                }
                self.absorb_internal(&mut cs, absorb_index, elements)?;
            },
            // We are in a duplex sponge construction. 
            // It is safe to squeeze the rate and absorb the new input without permuting 
            SpongeState::Squeezing{ next_squeeze_index } => { 
                self.absorb_internal(&mut cs, 0, elements)?;
            },
        };
        Ok(())
    }

    // Absorbs everything in elements, this does not end in an absorbtion.
    fn absorb_internal<CS: ConstraintSystem<F>>(&mut self, mut cs: CS, rate_start_index: usize, elements: &[FG]) -> Result<(), SynthesisError>
    {
        // if we can finish in this call
        if rate_start_index + elements.len() <= self.rate 
        {
            for i in 0..elements.len()
            {
                self.state[i + rate_start_index].add_in_place(&mut cs, &elements[i])?;
            }
            self.mode = SpongeState::Absorbing{ next_absorb_index: rate_start_index + elements.len() };
            return Ok(());
        }
        // otherwise absorb (rate - rate_start_index) elements
        let num_elements_absorbed = self.rate - rate_start_index;
        for i in 0..num_elements_absorbed
        {
            self.state[i + rate_start_index].add_in_place(&mut cs, &elements[i])?;
        }
        self.permute(&mut cs)?;
        // Tail recurse, with the input elements being truncated by num elements absorbed
        self.absorb_internal(cs, 0, &elements[num_elements_absorbed..])
    }

    pub fn squeeze<CS: ConstraintSystem<F>>(&mut self, mut cs: CS, num_elements: usize) -> Result<Vec<FG>, SynthesisError>
    {
        let mut squeezed_elems = AlgebraicSpongeGadget::<F, FG, P>::zero_state(&mut cs, num_elements)?;
        match self.mode {
            SpongeState::Absorbing{ next_absorb_index } => { 
                self.permute(&mut cs.ns(|| "absorb to squeeze"))?;
                self.squeeze_internal(&mut cs, 0, &mut squeezed_elems)?;
            },
            SpongeState::Squeezing{ next_squeeze_index } => {
                let mut squeeze_index = next_squeeze_index;
                if (squeeze_index == self.rate)
                {
                    self.permute(&mut cs)?;
                    squeeze_index = 0;
                }
                self.squeeze_internal(&mut cs, squeeze_index, &mut squeezed_elems)?;
            },
        };

        Ok(squeezed_elems)
    }

    // Squeeze |output| many elements. This does not end in a squeeze
    fn squeeze_internal<CS: ConstraintSystem<F>>(&mut self, mut cs : CS, rate_start_index: usize, output: &mut[FG]) -> Result<(), SynthesisError>
    {
        // if we can finish in this call
        if rate_start_index + output.len() <= self.rate 
        {
            for i in 0..output.len()
            {
                output[i] = self.state[i + rate_start_index].clone();
            }
            self.mode = SpongeState::Squeezing{ next_squeeze_index: rate_start_index + output.len() };
            return Ok(());
        }
        // otherwise squeeze (rate - rate_start_index) elements
        let num_elements_squeezed = self.rate - rate_start_index;
        for i in 0..num_elements_squeezed
        {
            output[i] = self.state[i + rate_start_index].clone();
        }

        // Unless we are done with squeezing in this call, permute.
        if (output.len() != self.rate)
        {
            self.permute(&mut cs)?;
        }
        // Tail recurse, with the correct change to indices in output happening due to changing the slice
        self.squeeze_internal(cs, 0, &mut output[num_elements_squeezed..])
    }
}
use algebra::prelude::*;
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::prelude::*;


#[derive(Copy, Clone)]
pub struct Domain<F>
where
    F: PrimeField,
{
    pub gen : F,
    pub offset : F,
    pub dim : u64,
}

impl<F : PrimeField> Domain<F> {

    pub fn order(&self) -> u64 
    {
        1 << self.dim
    }

    // Returns g, g^2, ... g^{dim}
    fn powers_of_gen(&self, dim : usize) -> Vec<F> {
        let mut result = Vec::new();
        let mut cur = self.gen;
        for i in 0..dim {
            result.push(cur);
            cur = cur * cur;
        }
        result
    }

    pub fn query_position_to_coset<CS: ConstraintSystem<F>, FG>(&self, 
        mut cs: CS,
        query_pos: &[Boolean],
        coset_dim: u64) -> Result<Vec<FG>, SynthesisError>
        where FG: FieldGadget<F, F>
    {
        let mut coset_index = query_pos;
        assert!(query_pos.len() == self.dim as usize || query_pos.len() == (self.dim - coset_dim) as usize);
        if query_pos.len() == self.dim as usize {
            coset_index = &coset_index[0..(coset_index.len() - coset_dim as usize)];
        }
        let mut coset = Vec::new();
        let powers_of_g = &self.powers_of_gen(self.dim as usize)[(coset_dim as usize)..];

        // build x = h * sum b_i * g^i
        let mut first_point_in_coset = FG::zero(&mut cs)?;
        // do summation
        let zero = FG::zero(&mut cs.ns(|| "Wasted zero allocation"))?;
        for i in 0..coset_index.len() {
            let term = zero.conditionally_add_constant(&mut cs.ns(|| format!("g^{:?} coefficient", i)), &coset_index[i], powers_of_g[i])?;
            first_point_in_coset.add_in_place(&mut cs, &term)?;
        }

        // multiply by domain offset
        first_point_in_coset.mul_by_constant_in_place(&mut cs, &self.offset);

        // build coset
        coset.push(first_point_in_coset);
        for i in 1..(1 << (coset_dim as usize))
        {
            let new_elem = coset[i - 1].mul_by_constant(&mut cs.ns(|| format!("coset term {:?}", i)), &self.gen)?;
            coset.push(new_elem);
        }

        Ok(coset)
    }

}
use algebra::fields::{PrimeField, Field};
use r1cs_core::{ConstraintSystem, SynthesisError};
use r1cs_std::{prelude::*};
use crate::hashing::{*, sponge::*, leafhash::*, two_to_one_hash::*};
use crate::algebra::mux::*;

// TODO: Make LeafHash and NodeHash traits
pub struct MerklePathGadget<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG> + Clone>
{
    pub leafHasher: LeafHashGadget<F, FG, P>,
    pub nodeHasher: TwoToOneHashGadget<F, FG, P>,
}

impl<F: PrimeField, FG: FieldGadget<F, F>, P: PermutationGadget<F, FG> + Clone> MerklePathGadget<F, FG, P>
{
    // We don't use the Zexe API for MTs.
    // The Zexe API for an authentication path requires the caller to pass in both the left and right hash
    // at every level. This is both computationally wasteful, and a pain for the caller to use.
    // Furthermore it requires a boolean should_enforce, which is a concern that can be handled 
    // at a higher level of abstraction.
    //
    // The path passed in should be the index of the leaf in the MT, represented in ?-endian bits. 
    // The trailing bit of the Path should be index % 2.
    // TODO: Recall which endianness is which
    // 
    // The authentication path should be given according to the same endianness.
    //  <TODO: Add proper explanation and diagram>, import prior diagram from Tendermint
    pub fn verify<CS: ConstraintSystem<F>>(&self, 
        mut cs: CS, 
        path: &[Boolean], 
        root: FG, 
        auth_path: &[FG], 
        leaf: &[FG]) -> Result<(), SynthesisError>
    {
        let leafHash = self.leafHasher.hash(&mut cs.ns(|| "leafhash"), leaf)?;
        let mut curHash = leafHash;

        // To traverse up a MT, we iterate over the path in reverse.
        let mut rev_path = path.to_vec();
        rev_path.reverse();
        let mut auth_path_rev = auth_path.to_vec();
        auth_path_rev.reverse();
        assert_eq!(auth_path_rev.len(), path.len());

        for i in (0..auth_path_rev.len()) {
            let mut cs_i = cs.ns(|| format!("MT layer {:?}", i));
            let bit = rev_path[i];
            // At any given bit, the bit being 0 indicates our currently hashed value is the left,
            // and the bit being 1 indicates our currently hashed value is on the right.
            // so we set left_hash = mux(bit, [curHash, authPathHash]), and similarly,
            // right_hash = mux(not bit, [curHash, authPathHash])
            // It is straight forward to generalize this to higher arity trees

            let curVec = vec![curHash.clone(), auth_path_rev[i].clone()];
            let left_hash = mux(&mut cs_i.ns(|| "left mux"), 
                &curVec, &[bit])?;
            let right_hash = mux(&mut cs_i.ns(|| "right mux"), 
                &curVec, &[bit.not()])?;

            curHash = self.nodeHasher.hash(&mut cs_i, left_hash, right_hash)?;
        }

        curHash.enforce_equal(&mut cs, &root)?;
        Ok(())
    }
}


#[cfg(test)]
mod test {
    use r1cs_std::{prelude::*, test_constraint_system::TestConstraintSystem};
    use r1cs_core::ConstraintSystem;
    use crate::merkle_tree::*;
    use crate::hashing::{two_to_one_hash::*, leafhash::*, dummy_permutation::*};
    use crate::algebra::{domain::Domain, mux::*, polynomial::DensePolynomial};

    use algebra::{FpParameters, prelude::*};
    use r1cs_std::prelude::*;
    use r1cs_core::SynthesisError;
    use r1cs_std::alloc::*;
    use r1cs_std::eq::EqGadget;
    use crate::alt_bn128::fr_gadget::FrGadget;
    use crate::alt_bn128::fr::Fr;
    use std::str::FromStr;

    #[test]
    fn mt_test() -> Result<(), SynthesisError> {
        let mut cs = TestConstraintSystem::<Fr>::new();
        let mut leaf_permutation = SeededDummyPermutation{seed: Fr::one()};
        let mut two_to_one_permutation = SeededDummyPermutation{seed: Fr::one() + Fr::one()};

        let leafsponge = AlgebraicSpongeGadget::<Fr, FrGadget, SeededDummyPermutation<Fr>>::
            new(&mut cs.ns(|| "1"), 1, 1, leaf_permutation)?;

        let nodesponge = AlgebraicSpongeGadget::<Fr, FrGadget, SeededDummyPermutation<Fr>>::
            new(cs.ns(|| "2"), 1, 1, two_to_one_permutation)?; 

        let leaf_hash = LeafHashGadget::new(&mut cs, leafsponge);
        let node_hash = TwoToOneHashGadget::new(&mut cs, nodesponge);

        // MT of {1, 2, 3, 4}
        let one = FrGadget::alloc(&mut cs.ns(|| format!("generate_{:?}", 1)), || Ok(Fr::from(1u32)))?;
        let two = FrGadget::alloc(&mut cs.ns(|| format!("generate_{:?}", 2)), || Ok(Fr::from(2u32)))?;
        let three = FrGadget::alloc(&mut cs.ns(|| format!("generate_{:?}", 3)), || Ok(Fr::from(3u32)))?;
        let four = FrGadget::alloc(&mut cs.ns(|| format!("generate_{:?}", 4)), || Ok(Fr::from(4u32)))?;
        // l1 = Leaf(1, 2) = 5 with the given permutation and seed
        let l1 = leaf_hash.hash(&mut cs.ns(|| "leaf hash 1"), &[one.clone(), two.clone()])?;
        assert_eq!(l1.get_value().unwrap(), Fr::from(5u32));
        // l2 = Leaf(3, 4) = 10
        let l2 = leaf_hash.hash(&mut cs.ns(|| "leaf hash 2"), &[three, four])?;
        assert_eq!(l2.get_value().unwrap(), Fr::from(9u32));
        // inner(l1, l2) = 22
        let root = node_hash.hash(&mut cs.ns(|| "node hash"), l1, l2.clone())?;

        let MT_path_gadget = MerklePathGadget{leafHasher: leaf_hash, nodeHasher: node_hash};
        MT_path_gadget.verify(&mut cs.ns(|| "MT 1"), &[Boolean::constant(false)], root, &[l2], &[one, two])?;

        assert!(cs.is_satisfied());
        Ok(())
    }
}
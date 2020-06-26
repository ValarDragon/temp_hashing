A copy of some of the optimizations we made in libiop, which are relevant to constraint complexity

* Optimal near-MDS matrix: 1/n percent of MDS entries are 0
* Alpha=17, lower the number of partial rounds significantly, and out-of-circuit eval time
* Poseidon paper, appendix A: Don't pay for ARKs or full-MDS in partial rounds. Push them all to the round before the partial rounds.

General note, we should be sure that we can balance the constraints evenly. Per Coda's work, the generic techniques for balancing non-zero terms across matrices allows to ~perfectly balance the Poseidon hash entries across all three matrices.

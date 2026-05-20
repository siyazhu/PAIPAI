# finiteT Mode Design

## Goal

finiteT mode is intended for proper finite-temperature Monte Carlo sampling.

The purpose is:
- Equilibrium sampling
- Thermodynamic averages
- Statistical properties
- Finite-temperature configurations

It is NOT intended for aggressive global optimization.

---

# Core Principles

finiteT mode should maintain:
- Markov-chain behavior
- Current-state-based acceptance
- Detailed-balance-like logic

Do NOT compare trial structures against:
- Global best structures
- Historical minimum energies
- Parallel-chain best states

---

# Acceptance Logic

Each move should compare:
- Proposed state
- Current state only

Typical Metropolis logic:
- Always accept if lower energy
- Otherwise probabilistic acceptance

---

# Cluster Interstitial Moves

finiteT mode may use `clusterInterstitialSwap` as an ordinary proposal move.

For thermodynamic correctness:
- generate one random cluster move at a time
- do not generate several local shuffles and choose the lowest-energy one
- accept or reject the relaxed proposal only against the current state

The current implementation swaps two different interstitial occupations and randomly shuffles the metal atom types in the union of their interstitial-site metal neighbors.

---

# Initial Structure

Recommended workflow:
1. Start from initial structure
2. Relax once using slow worker
3. Accept unconditionally
4. Begin normal MC sampling

---

# Relaxation Accuracy

finiteT mode may require:
- More accurate relaxation
- Smaller fmax
- Better geometry convergence

Loose relaxations may distort thermodynamic sampling.

---

# Drift Correction

Long simulations may exhibit:
- Global translation drift

Correction:
- Compute average metal displacement
- Remove overall translation

This preserves structural statistics.

---

# Worker Recommendations

finiteT mode generally prefers:
- Slow/refine workers only

Fast screening logic is discouraged.

GPU acceleration already provides sufficient performance.

---

# Important Risks

## Over-Optimization

Too much geometry optimization may:
- Artificially suppress entropy
- Distort finite-temperature behavior

Balance efficiency and physical correctness carefully.

---

## Interstitial Migration

Interstitials may migrate between cages during relaxation.

This may:
- Corrupt empty-site reconstruction
- Produce invalid site mappings

Special care is required.

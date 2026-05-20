# Search Mode Design

## Goal

Search mode is designed for:
- Efficient low-energy structure discovery
- Ground-state-like searches
- Broad configurational exploration

Thermodynamic correctness is NOT the primary objective.

---

# Main Strategy

Search mode uses:
- Fast workers
- Slow workers
- Multi-worker competition
- Screening logic

The purpose is to maximize exploration efficiency.

---

# Fast Workers

Purpose:
- Rapid approximate relaxation
- High throughput

Typical settings:
- FIRE optimizer
- Loose convergence

---

# Slow Workers

Purpose:
- Accurate refinement
- Final structure validation

Typical settings:
- LBFGS
- Tight convergence

---

# Cluster Interstitial Moves

Search mode may use the optional `clusterInterstitialSwap` move.

This move:
- swaps the occupations of two interstitial sites with different occupations
- collects the union of their neighboring metal atoms from `intsite_metal_neighbors.dat`
- randomly shuffles the metal atom types in that local neighbor union once

The purpose is to help interstitial atoms move together with a partially rearranged local metal cage, which can reduce trapping by strongly favorable local environments.

---

# Acceptance Philosophy

Search mode may:
- Prefer lower-energy structures aggressively
- Use global-best comparisons
- Prioritize exploration efficiency

Detailed balance is not required.

---

# Parallelism

Search mode is intended to scale well with:
- Multiple CPU cores
- Multiple GPUs
- Large worker pools

---

# Future Possibilities

Potential future improvements:
- Parallel tempering
- Block moves
- Adaptive moves
- Local rearrangements
- Exchange moves

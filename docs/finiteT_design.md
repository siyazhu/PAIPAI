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

If the initial relaxation reaches `--max-steps-refine` before satisfying `--fmax-refine`, PAIPAI still accepts it as the initial state as long as the worker returns a final structure and energy. The worker records `relax_converged` and `relax_max_force` in its metadata so the startup quality remains inspectable.

---

# Resume From Search

finiteT can start from a state previously accepted during search:

```bash
paipai --mode finiteT --resume-state search/mcprocess/000031 --temp 700
```

The resume directory should contain:
- `SAVE`: reference lattice and occupation after any search-mode reassignment
- `CONTCAR`: relaxed physical coordinates
- `energy` or `meta.json`: current relaxed energy

When `--resume-state` is used, PAIPAI skips the initial relaxation and begins finiteT sampling from that state. If `--root` is not specified, the launcher uses `finiteT_<temp>` as the working directory.

Use `SAVE`, not `REFERENCE_SAVE`, as the resumed occupation state. `REFERENCE_SAVE` records the pre-relaxation trial, while `SAVE` records the reconciled accepted state.

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

finiteT treats this as an invalid proposal for the discrete interstitial-site ensemble.

After relaxation, each interstitial atom is assigned to a candidate site using the cage-updated site positions. The assignment distance is limited by:

```bash
--interstitial-site-cutoff
```

If an interstitial atom is closest to a different site than the one occupied in the trial reference state, the proposal is rejected and the chain remains at the previous state.

The rejection is recorded in `mc.log` as:

```text
REJECT_INTERSTITIAL_HOP
```

This keeps the sampled state definition clear:

```text
state = reference interstitial-site occupation
energy = relaxed energy without changing site identity
```

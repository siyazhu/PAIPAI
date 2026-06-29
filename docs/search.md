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

# Step Budget

In search mode, `--steps N` is interpreted as a fast-screening budget.

The master stops once `counters/fast_count` reaches `N`, meaning about `N` structures have been processed by fast workers. This is different from `finiteT`, where `--steps` counts processed MC proposals.

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

# Local Hop Moves

Explicit `hop_interstitial` proposals are not recommended for search/prefast runs. During relaxation, these trials can be reassigned back to the original site pair or to a third site, which makes the trial proposal descriptor inconsistent with the final relaxed energy label. Use `--p-hop-inter 0` unless deliberately testing hop proposals.

---

# Prefast Warmup

When prefast is enabled, `--prefast-warmup-steps N` uses the first `N` valid MC proposals for learning without multi-candidate ranking. During warmup, each free fast-worker slot receives one ordinary random trial. The relaxed result still updates prefast weights using the final reassigned structure.

Discarded no-change tasks do not increment the warmup counter.

---

# Relaxed Interstitial Reassignment

Search mode allows interstitial atoms to hop between candidate sites during relaxation.

After slow refinement:
- relaxed interstitial atoms are matched to cage-updated candidate sites
- the assignment distance is limited by `--interstitial-site-cutoff`
- the shortest atom-site assignments are taken first
- one site can receive at most one interstitial atom

If the relaxed atom lands closer to a different valid site, PAIPAI updates the site occupation in the task `SAVE`. The task `CONTCAR` is then rewritten so its atom ordering matches the reassigned occupation state.

The reassignment is logged as:

```text
SEARCH_REASSIGN task_id=... n_reassigned=...
```

If no valid one-to-one assignment can be found within the cutoff, the candidate is discarded.

If reassignment leaves the final discrete occupation identical to the base state, PAIPAI discards the task without counting it as an MC proposal and without updating prefast. This is logged as:

```text
DISCARD_NO_CHANGE task_id=... move=... delta_source=... reason=... n_reassigned=...
```

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

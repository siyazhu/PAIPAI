# PAIPAI Architecture

## High-Level Overview

PAIPAI combines:
- Monte Carlo configurational sampling
- MLIP-based structural relaxation
- Persistent Python worker pools
- FIFO-based communication

The system consists of:
1. C++ Monte Carlo driver
2. Python relaxation workers
3. Structure storage and bookkeeping
4. Worker communication layer
5. Standalone input-building utilities

---

# Main Components

## C++ Driver

Responsibilities:
- Monte Carlo move proposals
- Acceptance/rejection decisions
- Worker scheduling
- Structure bookkeeping
- State management

---

## findinter Utility

Responsibilities:
- Read metal-only POSCAR files
- Generate or read effective hard-sphere radii from `radii.dat`
- Search for candidate interstitial sites on a periodic fractional grid
- Write PAIPAI `struc.in`
- Optionally write a POSCAR visualization with H markers at candidate sites

`findinter` is separate from the Monte Carlo driver. It prepares input structures before PAIPAI sampling begins.

---

## Python Workers

Responsibilities:
- Structure relaxation
- MLIP evaluation
- Force calculations
- Geometry optimization

Workers remain persistent throughout execution.

---

# Worker Communication

Communication uses:
- FIFO pipes
- Persistent workers

Advantages:
- Avoid repeated MLIP imports
- Avoid repeated Python startup overhead
- Better GPU utilization

---

# Worker Types

## Fast Workers

Purpose:
- Rapid screening
- Loose relaxation
- High throughput

Typical settings:
- FIRE optimizer
- Larger fmax

---

## Slow Workers

Purpose:
- Accurate refinement
- Final energy evaluation

If refinement reaches the configured maximum step count but still returns a final structure and energy, the worker writes that result normally and records `relax_converged` / `relax_max_force` in metadata. This lets the initial relaxation become the starting state even when it has not met the requested force tolerance.

Typical settings:
- LBFGS
- Smaller fmax

---

# Structure Flow

Typical workflow:

1. Generate MC trial structure
2. Send to worker
3. Relax structure
4. Return energy + coordinates
5. Reconcile relaxed interstitial atoms with reference sites
6. Accept/reject

---

# Reference And Relaxed States

PAIPAI keeps the discrete Monte Carlo state separate from the relaxed physical structure.

Files:
- `SAVE`: current reference lattice and site occupations
- `REFERENCE_SAVE`: reference state used to launch a specific relaxation task
- `CONTCAR`: relaxed physical coordinates from the worker

The reference site ordering in `SAVE` is fixed. `CONTCAR` is ordered by atom species, so PAIPAI maps relaxed atoms back to reference site ids using the same ordering convention used by `Structure::outputvasp()`.

When generating a new trial `POSCAR`, PAIPAI seeds coordinates from the accepted `CONTCAR` for efficiency:
- relaxed metal coordinates are reused
- already occupied interstitial sites reuse their relaxed interstitial coordinates when the site remains occupied
- newly occupied interstitial sites are placed by updating the target site from its relaxed metal cage

The trial `SAVE` still stores only the reference occupation state.

---

# Resume State

A finiteT run can resume from an accepted state directory, typically one under `search/mcprocess/`.

Required state files:
- `SAVE`
- `CONTCAR`
- `energy` or `meta.json`

The resumed `SAVE` is the accepted reference state. This is important for search outputs because `SAVE` may include relaxed interstitial reassignment, while `REFERENCE_SAVE` only records the trial before relaxation.

If `--resume-state` points at an `mcprocess` directory rather than one numbered state, PAIPAI selects the latest numbered child directory. Resume mode does not require an input structure file.

---

# Interstitial Site Reconciliation

After slow relaxation, relaxed interstitial atoms are assigned back to candidate sites by distance to the cage-updated site positions. The maximum assignment distance is controlled by:

```bash
--interstitial-site-cutoff
```

If multiple assignments are possible, PAIPAI greedily assigns the shortest atom-site distances first while enforcing one interstitial per site.

Mode-dependent behavior:
- `finiteT`: if any relaxed interstitial maps to a different site than the one in the trial reference state, the trial is rejected.
- `search`: relaxed hopping is allowed; the reference occupation is rewritten to the reassigned sites and `CONTCAR` is reordered to remain consistent with the new `SAVE`.

---

# Neighbor Map

Neighbor maps:
- Built once initially
- Stored in intsite_metal_neighbors.dat

Purpose:
- Track local cage deformation
- Move empty sites consistently

---

# finiteT Workflow

finiteT mode differs from search mode:
- Single Markov chain
- Strict current-state comparison
- Thermodynamic correctness prioritized
- Interstitial site hopping during relaxation is rejected

---

# Search Workflow

Search mode:
- Parallel structure exploration
- Global-best competition
- Efficiency-focused
- Relaxed interstitial site hopping is accepted by reassigning site occupations

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
5. Update SAVE
6. Accept/reject

---

# SAVE Update Logic

After relaxation:
- Metals updated directly
- Occupied interstitials follow CONTCAR
- Empty sites reconstructed from neighbor displacement

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

---

# Search Workflow

Search mode:
- Parallel structure exploration
- Global-best competition
- Efficiency-focused

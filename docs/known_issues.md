# Known Issues

## Interstitial Cage Migration

Problem:
- Interstitial atoms may migrate between cages during relaxation.

Potential consequence:
- Empty-site reconstruction may become invalid.
- Some interstitial sites may collapse unnaturally close together.

---

# Drift During finiteT Sampling

Problem:
- Structures may slowly drift spatially during long simulations.

Cause:
- Coordinate reconstruction accumulates translation errors.

Mitigation:
- Remove average metal displacement after relaxation.

---

# Overly Loose Relaxation

Problem:
- Large fmax values may produce inaccurate energies.

Potential consequence:
- Incorrect Monte Carlo statistics
- Distorted acceptance probabilities

---

# refine_outbox Growth

Problem:
- refine_outbox may accumulate large numbers of temporary files.

Mitigation:
- Clean processed refinement files automatically.

---

# GPU Resource Contention

Problem:
- Running many MLIP relaxations on the same GPU may reduce efficiency.

Mitigation:
- Limit concurrent heavy worker usage per GPU.

---

# Structure Mapping Assumptions

Current logic assumes:
- Structures originate from PAIPAI-generated ordering conventions.

External structures may break assumptions.

---

# Very Open Structures

Neighbor-based empty-site reconstruction may become unreliable for:
- Highly distorted systems
- Very open structures
- Large local rearrangements

---

# findinter Effective Radii

Problem:
- `findinter` uses effective hard-sphere radii, and the site count can change strongly with `radii.dat`.

Potential consequence:
- Too-large interstitial radii may miss valid tetrahedral or octahedral sites.
- Too-small radii may introduce extra candidate sites.

Mitigation:
- Inspect the optional H-marker POSCAR from `--site-poscar`.
- For periodic reference systems, compare the number of detected sites against the expected crystallographic count.

---

# findinter Grid Resolution

Problem:
- Very coarse grids may miss local maxima of the clearance field.

Mitigation:
- Use grids compatible with supercell repeats when possible.
- For example, a 5x5x5 BCC-like supercell is naturally tested with `--grid 100`.

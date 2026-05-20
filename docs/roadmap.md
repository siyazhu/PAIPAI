# PAIPAI Roadmap

## Short-Term Goals

- Improve finiteT robustness
- Improve coordinate reconstruction stability
- Better logging and diagnostics
- Cleaner file management
- Better GPU scheduling
- Further validate `findinter` on distorted defects, surfaces, and grain boundaries

---

# Possible Future Features

## Parallel Tempering

Purpose:
- Improve configurational sampling
- Escape local minima

---

## Block Moves

Purpose:
- Improve sampling efficiency
- Correlated configurational updates

---

## Adaptive Moves

Purpose:
- Dynamically adjust move probabilities
- Improve acceptance behavior

---

## Local Rearrangement

Purpose:
- Relax local environments after exchanges
- Improve physically meaningful moves

---

## Interstitial Site Detection

Purpose:
- Improve and validate the new `findinter` utility
- Improve generality

Current status:
- Initial grid-based `findinter` implementation is available
- Supports `radii.dat`, `min/max-void-factor`, PAIPAI `struc.in` output, and H-marker POSCAR visualization

Future work:
- Better defaults for effective radii
- Optional Voronoi-assisted candidate generation
- Optional large-void packing mode for storage-material workflows

---

## MD-Assisted Sampling

Purpose:
- Occasionally introduce MD-based perturbations
- Improve exploration

---

# Long-Term Goals

- Better thermodynamic integration
- Advanced defect sampling
- Improved free-energy calculations
- More generalized lattice support
- Multi-component interstitial systems
- Better distributed parallelism

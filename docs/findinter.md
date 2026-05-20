# findinter Design and Usage

`findinter` is a standalone PAIPAI utility for building `struc.in` files from a metal-only POSCAR.

It is intended for systems where interstitial sites are not already known, including:
- multi-component bulk alloys
- defective crystals
- grain boundaries
- surfaces
- other distorted metallic structures

---

# Basic Workflow

Typical usage:

```bash
findinter \
  --input POSCAR_metal_only \
  --inter B \
  --internum 4 \
  --grid 100 \
  --min-void-factor 0.9 \
  --max-void-factor 1.2 \
  --output struc.in \
  --site-poscar interstitial_sites.vasp
```

Inputs:
- `--input`: metal-only POSCAR
- `--inter`: comma-separated interstitial element list
- `--internum`: comma-separated number of occupied interstitial atoms

Outputs:
- PAIPAI `struc.in`
- optional visualization POSCAR with every candidate interstitial site marked as H

---

# radii.dat

`findinter` uses effective hard-sphere radii.

If `radii.dat` is not present, `findinter` creates one from built-in defaults and exits. The user should inspect or edit the file, then rerun the same command.

The radii are search parameters, not immutable tabulated atomic radii. They control how permissive the interstitial-site search is.

Example:

```text
Nb 1.428
Ta 1.428
Ti 1.428
Hf 1.428
Interstitial 0.55
```

For a 5x5x5 BCC-like NbTaTiHf cell with 250 metal atoms, this setup with `--grid 100 --min-void-factor 0.9 --max-void-factor 1.2` identifies 1500 tetrahedral candidate sites.

---

# Algorithm

`findinter` works in fractional grid space.

For each grid point:
1. Convert the fractional coordinate to Cartesian coordinates.
2. Compute the nearest metal distance using periodic boundary conditions.
3. Keep the point only if:

```text
d_nearest >= min_void_factor * (r_metal + r_interstitial)
d_nearest <= max_void_factor * (r_metal + r_interstitial)
```

The program then searches the clearance field for local maxima.

After selecting one interstitial site, all grid points within:

```text
2 * r_interstitial * min_void_factor
```

are suppressed to a very large negative clearance. The local maxima are then searched again. This prevents nearly overlapping candidate sites while still allowing nearby points to become local maxima after the closest selected site is removed.

---

# Choosing Grid Size

For periodic crystals and supercells, choose a grid compatible with the supercell repeat whenever possible.

For example, a 5x5x5 supercell works naturally with:

```text
--grid 100
```

because each conventional-cell repeat receives the same number of grid intervals.

For distorted, defective, or non-cubic cells, larger grids may improve site recovery but increase runtime roughly with `grid^3`.

---

# Important Limitations

`findinter` is a practical site generator, not a full Voronoi or global packing solver.

Known limitations:
- Effective radii strongly affect the resulting site count.
- Very coarse grids can miss narrow local maxima.
- The iterative hard-sphere suppression step can change the number of retained sites.
- Large open voids and channels are not globally packed or optimized.

For storage-material style workflows where the goal is to pack many H-like atoms into large voids, a future dedicated packing mode may be more appropriate.


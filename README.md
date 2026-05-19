# Package for Alloy Interstitial Predictions using Artificial Intelligence (PAIPAI)

<p>
PAIPAI is a versatile computational tool designed to efficiently search for crystalline metallic structures—especially those containing interstitials, point defects, grain boundaries, or surface slabs—with the lowest free energy.
By combining Monte Carlo sampling techniques with machine-learning interatomic potentials (MLIPs), PAIPAI enables rapid and accurate exploration of vast configurational spaces that are traditionally inaccessible to first-principles methods alone.
</p>

<p>
PAIPAI is particularly designed for:
</p>

- High-entropy alloys (HEAs)
- Interstitial solutes (O, B, C, N, etc.)
- Defect-containing structures
- Large-scale chemically disordered systems
- Ground-state structure search
- Finite-temperature Monte Carlo sampling

<p>
  <a href="https://github.com/siyazhu/PAIPAI/issues/new?labels=bug">Report a Bug</a> |
  <a href="https://github.com/siyazhu/PAIPAI/issues/new?labels=enhancement">Request a Feature</a>
</p>

---

# Versions

## Stable version (`main`)

The `main` branch contains the original stable PAIPAI implementation.

Recommended for:
- reproducing older calculations
- conservative production workflows
- users who prefer a more established version

---

## Development version (`v2.0-dev`)

The `v2.0-dev` branch contains the new PAIPAI v2.0 framework.

Major updates include:
- GPU-supported MLIP calculations
- new `search` and `finiteT` modes
- improved structure-update workflow after relaxation
- reduced relaxation cost
- cleaner output handling
- multiple bug fixes

The v2.0 branch is currently under active development and testing.

---

# Features

- MLIP-based structure optimizations and energy evaluations
- Efficient Monte Carlo simulations
- Parallel worker-based structure search
- Finite-temperature Markov-chain Monte Carlo sampling
- Supports complex multi-element metallic systems with interstitials and various defect-containing structures
- GPU-supported MLIP calculations
- Dual-worker search workflow for efficient ground-state exploration

---

# Repository Structure

```text
PAIPAI/
├── examples/
├── include/
├── scripts/
├── src/
├── CMakeLists.txt
└── README.md
```

- `src/` : C++ source files
- `include/` : C++ header files
- `scripts/` : Python workers and PAIPAI launch scripts
- `examples/` : example input files and SLURM scripts

---

# Installation

## Dependencies

PAIPAI relies on the external MLIP infrastructure provided by MaterialsFramework:

https://github.com/dogusariturk/MaterialsFramework

Please install MaterialsFramework and ensure the corresponding Python environment is properly configured before using PAIPAI.

Typical dependencies include:
- pymatgen
- TensorFlow
- MLIP models
- materialsframework package

Users should verify that the following command works correctly before running PAIPAI:

```bash
python -c "import materialsframework"
```
---

## Prerequisites

- CMake ≥ 3.10
- C++17 compatible compiler (e.g. GCC, Clang)
- Python 3
- MaterialsFramework

---

## Clone repository

```bash
git clone https://github.com/siyazhu/PAIPAI.git
cd PAIPAI
```

To use the development version:

```bash
git checkout v2.0-dev
```

---

# Build and Install (PAIPAI v2.0)

```bash
mkdir build
cd build

cmake .. -DCMAKE_INSTALL_PREFIX=$HOME

cmake --build .

cmake --install .
```

This installs:

```text
$HOME/bin/paipai
$HOME/bin/mc_paipai
$HOME/bin/fast_worker.py
$HOME/bin/slow_worker.py
```

If necessary, add `$HOME/bin` to your PATH:

```bash
export PATH="$HOME/bin:$PATH"
```

You may add this line to your `~/.bashrc`.

---

# Legacy Compilation (Original Version)

The original implementation can still be compiled manually:

```bash
g++ -std=c++17 -O2 -pthread \
constants.cpp mcpaipai.cpp element.cpp structure.cpp \
-o mc_paipai -lstdc++fs
```

---

# Basic Usage

After installation:

```bash
paipai --help
```

---

# Running PAIPAI

## Search mode

`search` mode uses a dual-worker architecture:
- fast workers perform coarse screening relaxations
- slow workers perform refinement relaxations

This mode supports many parallel workers and is mainly intended for:
- ground-state search
- low-temperature structure optimization
- rapid screening

Typical usage:

```bash
paipai \
  --input struc.in \
  --mode search \
  --device cpu \
  --fast 15 \
  --slow 15 \
  --steps 2000 \
  --temp 10
```

A relatively low Monte Carlo temperature is generally recommended for structure search.

---

## Finite-temperature mode

`finiteT` mode performs proper sequential Markov-chain Monte Carlo sampling.

Characteristics:
- uses only one slow worker
- does NOT support parallel MC chains
- generates unbiased finite-temperature ensembles

Typical usage:

```bash
paipai \
  --input struc.in \
  --mode finiteT \
  --device cuda \
  --ngpu 1 \
  --fast 0 \
  --slow 1 \
  --steps 200000 \
  --temp 700
```

---

# GPU Support

PAIPAI v2.0 supports GPU-accelerated MLIP relaxation and energy evaluation.

To use GPU acceleration:

```bash
--device cuda
```

Users may refer to the scripts in:

```text
examples/
```

for example GPU submission scripts.

---

# Command-Line Options

## Basic options

| Option | Description |
|---|---|
| `--input FILE` | Input structure file |
| `--mode MODE` | `search` or `finiteT` |
| `--device DEV` | `cpu` or `cuda` |
| `--model NAME` | MLIP model |
| `--root DIR` | Working directory |

---

## Worker options

| Option | Description |
|---|---|
| `--fast N` | Number of fast workers |
| `--slow N` | Number of slow workers |
| `--ngpu N` | Number of GPUs used for round-robin worker assignment |
| `--pool-cap N` | Waiting-pool capacity |

---

## Monte Carlo options

| Option | Description |
|---|---|
| `--steps N` | Number of Monte Carlo steps |
| `--temp T` | Monte Carlo temperature |
| `--p-swap-metal N` | Weight for metal swap moves |
| `--p-swap-inter N` | Weight for interstitial swap moves |
| `--intsite-neighbor-cutoff X` | Cutoff for interstitial-site neighbor mapping |

---

# Major Updates in v2.0

## 1. MLIP-based structure-update workflow

PAIPAI v2.0 introduces a new workflow for updating structural information after MLIP relaxation.

### (1) Metal-coordinate update

After each relaxation:
- relaxed metal coordinates are written back into the structure file
- subsequent relaxations therefore start closer to equilibrium structures

This significantly reduces the average number of optimization steps.

### (2) Interstitial-site coordinate update

For interstitial sites:
- neighboring metal atoms are determined using a cutoff distance
- the neighbor mapping is stored in:

```text
intsite_metal_neighbors.dat
```

After each relaxation:
- updated metal coordinates are used to reconstruct interstitial-site coordinates

This greatly improves relaxation efficiency and structural continuity.

However, for relatively open or low-density structures, additional testing is still needed to fully validate robustness.

---

## 2. New `search` and `finiteT` modes

PAIPAI v2.0 now supports two distinct workflows.

### `search` mode
- dual-worker architecture
- many-worker parallel execution
- optimized for structure search and ground-state exploration

### `finiteT` mode
- sequential Markov-chain Monte Carlo
- single slow-worker workflow
- unbiased finite-temperature ensemble sampling

---

## 3. GPU support

PAIPAI v2.0 supports:
- GPU-accelerated MLIP structural relaxation
- GPU-accelerated energy evaluation

GPU usage can be enabled using:

```bash
--device cuda
```

Worker processes can also be distributed across multiple GPUs using round-robin assignment.

---

## 4. Bug fixes and workflow cleanup

PAIPAI v2.0 additionally:
- fixes overall atomic drift issues
- simplifies output files
- cleans `refine_outbox` and temporary outputs more aggressively
- improves worker scheduling
- fixes several miscellaneous bugs

---

# Examples

Example scripts and input files are provided in:

```text
examples/
```

including:
- CPU search examples
- GPU finite-temperature examples
- SLURM submission scripts
- Sample struc.in file:
    - struc\_TiVCrRe\_Slab.in: input structure for Ti19V77Cr26Re6 slab structure;
    - struc\_NbTaTiHf\_Bulk\_B.in: input structure for Nb113Ti63Ta37Hf37 bulk structure, with 8 Boron interstitials;
    - struc\_NbTaTiHf\_GB\_BO.in: input structure for Nb90Ti50Ta30Hf30 sigma-5(120) grain boundary structure, with 4 Boron and 4 Oxygen interstitials. 
---

# Citation

If you use PAIPAI in academic work, please cite the associated publications: Zhu, Siya, and Raymundo Arróyave. "Ground-state structure search of defective high-entropy alloys using machine-learning potentials and Monte Carlo sampling." Computational Materials Science 270 (2026): 114752. https://doi.org/10.1016/j.commatsci.2026.114752

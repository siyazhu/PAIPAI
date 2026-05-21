<img width="200" alt="logo" src="logo.png"/>

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

PAIPAI supports multi-stage workflows including:

- Random or symmetry-constrained Monte Carlo sampling
- Machine-learning potential based structural relaxation
- Low-energy candidate filtering
- Free-energy-aware finite-temperature exploration
- Parallel execution on HPC clusters

![Workflow diagram](docs/workflow.png)

## Key Features

- **Flexible structure input** through CIF / POSCAR / ASE-readable files
- **Monte Carlo search** for substitutional and interstitial configurations
- **Support for multiple ML potentials**, including MACE, CHGNet, M3GNet, SevenNet, ORB, UMA, and others
- **Parallel task scheduling** for high-throughput local relaxation
- **Two-stage relaxation workflows** for efficient pre-screening and refinement
- **Finite-temperature Monte Carlo mode** with configurable Boltzmann acceptance
- **Detailed reports** for energies, accepted moves, and candidate structures

## Installation

Clone the repository:

```bash
git clone https://github.com/siyazhu/PAIPAI.git
cd PAIPAI
```

Create a build directory:

```bash
mkdir build
cd build
```

Configure with CMake:

```bash
cmake ..
```

Build:

```bash
make -j
```

Install:

```bash
make install
```

The main executable script will be installed as:

```bash
paipai
```

## Python Environment

PAIPAI relies on Python worker scripts for ML-potential-based relaxations. The Python environment should include packages such as:

- `ase`
- `numpy`
- `pymatgen`
- The selected ML potential backend, for example `mace`, `chgnet`, `m3gnet`, `sevenn`, `orb`, or `uma`

A typical conda environment may look like:

```bash
conda create -n paipai python=3.10
conda activate paipai
pip install ase numpy pymatgen
```

Then install your desired ML potential package following its official instructions.

## Basic Usage

Run PAIPAI with:

```bash
paipai --input paipai.in
```

A typical input file defines:

- Structure file
- Sampling mode
- Element types
- Number of Monte Carlo steps
- Temperature, if using finite-temperature mode
- ML potential backend
- Relaxation settings
- Output directory

Example:

```text
MODE search
STRUCTURE input.cif
ELEMENTS Nb O
N_STEPS 1000
TEMPERATURE 1000
CALCULATOR mace
OUTPUT_DIR run_paipai
```

## Modes

PAIPAI currently supports two main workflows:

### Search Mode

Search mode is designed to rapidly explore large configuration spaces and identify low-energy structures.

```text
MODE search
```

In this mode, PAIPAI proposes trial configurations, relaxes them using ML potentials, and stores low-energy candidates for further analysis.

### Finite-Temperature Mode

Finite-temperature mode performs Monte Carlo sampling using a Boltzmann-like acceptance rule.

```text
MODE finiteT
TEMPERATURE 1000
```

This is useful for studying thermally accessible configurations and approximate finite-temperature chemical disorder.

## Output

A typical PAIPAI run creates output directories such as:

```text
fast/
waiting_pool/
waiting_work/
refine_outbox/
reports/
```

Depending on the selected mode, PAIPAI records relaxed structures, energy summaries, accepted moves, and candidate configurations.

## Documentation

Additional documentation and examples will be added under:

```text
docs/
examples/
```

## Citation

If you use PAIPAI in your work, please cite the relevant PAIPAI publication or repository.

## License

Please see the repository license file for details.

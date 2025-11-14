# Package of Alloy Interstitials Processed with Artificial Intelligence (PAIPAI)


<p>
PAIPAI is a versatile computational tool designed to efficiently search for crystalline metallic structures—especially those containing interstitials, point defects, grain boundaries, or surface slabs—with the lowest free energy.
By combining Monte Carlo sampling techniques with machine-learning interatomic potentials (MLIPs), PAIPAI enables rapid and accurate exploration of vast configurational spaces that are traditionally inaccessible to first-principles methods alone.</p>

<p>
  <a href="https://github.com/siyazhu/paipai/issues/new?labels=bug">Report a Bug</a> |
  <a href="https://github.com/siyazhu/paipai/issues/new?labels=enhancement">Request a Feature</a>
</p>

</div>

---

## Features

- MLIP-based structure optimizations and energy evaluations
- Efficient Monte Carlo simulations
- Supports complex multi-element metallic systems with interstitials and various defect-containing structures.
---

## Installation

### Prerequisites

- CMake ≥ 3.10
- C++17 compatible compiler (e.g., GCC, Clang)
- Python 3 

### Build Instructions

```sh
git clone https://github.com/siyazhu/PAIPAI.git
cd PAIPAI
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME
cmake --build . --target install
```

---

## Quick Start

Use the `paipai` command with an input structure file (default: struc.in) to start a PAIPAI run.

```sh
paipai --input [Input Struc] --model [Model Name] --device [cpu/cuda] --fast [# fast workers] --slow [# slow workers] --steps [# MC steps] --temp [MC temperature] --pool-cap [Waiting pool capacity] 
```
or
```sh
paipai -h
```
for help information.
A sample struc.in file is in the /examples/.

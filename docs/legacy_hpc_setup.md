# Legacy HPC deployment notes (CentOS 7 / PBS)

These notes summarize a working PAIPAI + MaterialsFramework + GRACE deployment tested on an older CPU-only HPC system with:

- CentOS 7 / Linux 3.10
- glibc 2.17
- Bash 4.2.46
- system Python 2.7.5
- system GCC/G++ 4.8.5
- Torque/PBS (`qsub`, `qstat`)
- CPU nodes exposed as `q0`, typically `np = 32`

## Python and MaterialsFramework

Do not modify the system Python. A user-space Python 3.12 environment is required by current MaterialsFramework. `uv` can be used to create/manage the environment.

On this old platform, several packages may build from source because modern wheels are unavailable for glibc 2.17. The deployment required a newer user-space GCC/G++ plus development libraries including JPEG, OpenBLAS, CMake/pkg-config, and HDF5.

A working layout was:

```text
~/miniconda3/                  # bootstrap + modern compiler/runtime libraries
~/envs/materialsframework/    # Python 3.12 venv
```

The system `libstdc++.so.6` is too old for the resulting NumPy/TensorFlow binaries. The newer Conda runtime must be loaded for PAIPAI/GRACE jobs, e.g.

```bash
export LD_PRELOAD="$HOME/miniconda3/lib/libstdc++.so.6"
```

Do not replace `/lib64/libstdc++.so.6` system-wide.

## GRACE model cache on restricted networks

TensorPotential normally downloads `GRACE-2L-OMAT` from Hugging Face. If the HPC can access GitHub but cannot connect to Hugging Face, download the model elsewhere and place the extracted model at:

```text
~/.cache/grace/GRACE-2L-OMAT/
```

The model directory should contain the SavedModel files (`saved_model.pb`, `variables/`, `assets/`, metadata, etc.). A manually populated cache may trigger a provenance/update warning. If the model provenance has been independently verified, the update check can be suppressed with:

```bash
export GRACE_NO_UPDATE_CHECK=1
```

## MaterialsFramework 1.x compatibility

Current `GraceCalculator` no longer accepts PAIPAI's legacy `device=` keyword. Passing it through `**kwargs` eventually reaches `object.__init__()` and fails with:

```text
TypeError: object.__init__() takes exactly one argument (the instance to initialize)
```

Device selection should therefore not be passed as `GraceCalculator(device=...)`. For CPU-only batch jobs, device visibility can be controlled externally, for example:

```bash
export CUDA_VISIBLE_DEVICES=""
```

Current MaterialsFramework also uses `steps` for the relaxation step limit rather than the older `max_steps` keyword/property.

## Bash 4.2 + `set -u` empty-array compatibility

On Bash 4.2, an initialized-but-empty array can still produce an `unbound variable` failure when expanded under `set -u`, for example when optional master arguments are empty.

The robust workaround is to disable nounset locally in the subshell that launches the C++ master:

```bash
(
    cd "$ROOT_DIR"
    set +u
    "$MASTER" "${MC_EXTRA_ARGS[@]}" ...
)
```

Because this is a subshell, the outer PAIPAI wrapper retains `set -u` after the master invocation exits.

## PBS example (CPU-only)

A 32-slot q0 node can be used as 30 PAIPAI workers plus one master slot. Example:

```bash
#!/usr/bin/env bash
#PBS -N PAIPAI_search_cpu
#PBS -q q0
#PBS -l nodes=1:ppn=31
#PBS -l walltime=72:00:00
#PBS -j oe

set -eo pipefail
cd "$PBS_O_WORKDIR"

source "$HOME/envs/materialsframework/bin/activate"
export LD_PRELOAD="$HOME/miniconda3/lib/libstdc++.so.6"
export PAIPAI_PYTHON="$HOME/envs/materialsframework/bin/python"
export PATH="$HOME/envs/materialsframework/bin:$PATH"

export CUDA_VISIBLE_DEVICES=""
export GRACE_NO_UPDATE_CHECK=1

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export TF_NUM_INTRAOP_THREADS=1
export TF_NUM_INTEROP_THREADS=1

paipai \
  --input struc.in \
  --mode search \
  --device cpu \
  --fast 15 \
  --slow 15 \
  --pool-cap 128 \
  --steps 2000 \
  --temp 700 \
  --p-swap-metal 70 \
  --p-swap-inter 30 \
  --p-cluster-inter 0 \
  --p-exch-metal 0 \
  --p-exch-inter 0 \
  --intsite-neighbor-cutoff 3.5 \
  > paipai_search_cpu.log 2>&1
```

The explicit redirection is useful on Torque/PBS systems where scheduler stdout/stderr may remain in a compute-node spool until the job finishes. Worker logs written directly by PAIPAI remain visible during execution.

## PBS commands

```bash
qsub run_paipai.pbs       # submit
qstat                      # show own jobs
qdel JOBID.master          # cancel
qstat -Q                   # queue summary
pbsnodes -a                # node/resource information
```

## Notes from the test

- TensorFlow may print a CUDA initialization error even in a CPU-only job, then successfully initialize XLA for `Host`. If execution continues with `StreamExecutor device (0): Host` and XLA compilation succeeds, the CPU backend is working.
- Slow workers atomically claim waiting-pool tasks via `os.rename`, so multiple slow workers can operate concurrently when multiple candidates are available.
- Avoid globally exporting the newer `LD_PRELOAD` in `.bashrc`; keep it scoped to the MaterialsFramework/PAIPAI environment so unrelated legacy software continues to use the system runtime.

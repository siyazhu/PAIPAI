#!/usr/bin/env python3
"""
Fast worker: screen-only relaxation into waiting_pool.
- Inputs (from C++ master via 'fast/' directory):
    fast/POSCAR{k}, fast/SAVE{k}, fast/.go_{k}
- Outputs:
    fast/.done_{k}   # just a 'done' signal so C++ knows it can reuse this slot
- It does NOT write any RESULT back to C++.
- It runs a coarse relaxation (screen) and pushes the result into:
    waiting_pool/<task_id>/{POSCAR, SAVE, meta.json}
- The waiting_pool has a soft capacity:
    - If the pool is full, the worker evicts the highest-energy entry
      ONLY if its own energy is lower; otherwise it drops this candidate.
- Counters:
    ROOT/counters/fast_count will be incremented for each job handled.
"""
print("FAST WORKER START", flush=True)

import os, sys, json, time, tempfile, signal, uuid
from pathlib import Path
from typing import Optional
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar
from materialsframework.calculators import GraceCalculator

# ---------- atomic helpers ----------
def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp_", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as f:
            f.write(text)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, str(path))
    finally:
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass

def atomic_write_json(path: Path, obj) -> None:
    atomic_write_text(path, json.dumps(obj, ensure_ascii=False, indent=2) + "\n")

def atomic_write_poscar(structure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp_", dir=str(path.parent))
    os.close(fd)
    Poscar(structure).write_file(tmp)
    with open(tmp, "rb") as rf:
        os.fsync(rf.fileno())
    os.replace(tmp, str(path))
    try:
        os.remove(tmp)
    except FileNotFoundError:
        pass

def atomic_copy(src: Path, dst: Path) -> None:
    """Copy a small file atomically (read all, then replace)."""
    if not src.exists():
        return
    dst.parent.mkdir(parents=True, exist_ok=True)
    data = src.read_bytes()
    fd, tmp = tempfile.mkstemp(prefix=".tmp_", dir=str(dst.parent))
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, str(dst))
    finally:
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass

def atomic_increment_counter(root: Path, name: str, retries: int = 50, sleep_s: float = 0.01):
    """Naive lock-based atomic integer counter."""
    cdir = root / "counters"
    cdir.mkdir(parents=True, exist_ok=True)
    cpath = cdir / name
    lock = cdir / (name + ".lock")
    for _ in range(retries):
        try:
            fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
            os.close(fd)
        except FileExistsError:
            time.sleep(sleep_s)
            continue
        try:
            old = 0
            if cpath.exists():
                try:
                    old = int(cpath.read_text().strip() or "0")
                except Exception:
                    old = 0
            atomic_write_text(cpath, f"{old+1}\n")
        finally:
            try:
                lock.unlink()
            except FileNotFoundError:
                pass
        return

# ---------- waiting_pool helpers ----------
def list_pool_with_energy(pool: Path):
    """Return list of (energy_screen, path) for tasks in waiting_pool."""
    out = []
    if not pool.exists():
        return out
    for d in pool.iterdir():
        if not d.is_dir() or d.name.startswith(".tmp_"):
            continue
        meta = d / "meta.json"
        try:
            E = float(json.loads(meta.read_text())["energy_screen"])
        except Exception:
            E = 1e99
        out.append((E, d))
    return out

def maybe_evict_for_capacity(pool: Path, capacity: int, my_E: float) -> bool:
    """
    If pool size >= capacity, evict the highest-energy task IFF my_E < E_max.
    Returns True if we can insert, False if we should drop this candidate.
    """
    cands = list_pool_with_energy(pool)
    if len(cands) < capacity:
        return True
    cands.sort(key=lambda x: x[0])  # ascending
    E_max, p_max = cands[-1]
    if my_E < E_max:
        try:
            for p in p_max.iterdir():
                try:
                    p.unlink()
                except IsADirectoryError:
                    pass
            p_max.rmdir()
            return True
        except Exception:
            # If removal fails (race), re-check size once
            cands2 = list_pool_with_energy(pool)
            return len(cands2) < capacity
    else:
        return False

print("FAST WORKER PYTHON MODULE LOADED", flush=True)

# ---------- main ----------
def main():
    print("FAST WORKER STARTED", flush=True)
    import argparse
    ap = argparse.ArgumentParser(description="Fast worker (screen-only into waiting_pool)")
    ap.add_argument("--slot", type=int, required=True, help="worker slot index")
    ap.add_argument("--root", type=str, default=".", help="root directory")
    ap.add_argument("--model", type=str, default="GRACE-2L-OMAT", help="model name")
    ap.add_argument("--device", type=str, default="cpu", choices=["cpu", "cuda"])
    ap.add_argument("--fmax_screen", type=float, default=0.10)
    ap.add_argument("--max_steps_screen", type=int, default=30)
    ap.add_argument("--pool_cap", type=int, default=128, help="waiting_pool capacity")
    args = ap.parse_args()

    k = args.slot
    ROOT = Path(args.root).resolve()
    FAST = ROOT / "fast"
    POOL = ROOT / "waiting_pool"

    # handshake files with C++ master
    poscar = FAST / f"POSCAR{k}"
    savef  = FAST / f"SAVE{k}"
    gof    = FAST / f".go_{k}"
    donef  = FAST / f".done_{k}"

    calc = GraceCalculator(
        model=args.model,
        device=args.device,
        verbose=True,
        fmax=args.fmax_screen,
        max_steps=args.max_steps_screen,
    )
    sys.stderr.write(f"[fast {k}] initialized @ {ROOT}\n")

    stop = {"flag": False}
    def _sigterm(*_): stop["flag"] = True
    signal.signal(signal.SIGINT, _sigterm)
    signal.signal(signal.SIGTERM, _sigterm)

    while not stop["flag"]:
        if not gof.exists():
            time.sleep(0.05)
            continue

        try:
            atomic_increment_counter(ROOT, "fast_count")

            # read POSCAR (retry a short window)
            struct = None
            for _ in range(10):
                try:
                    struct = Structure.from_file(str(poscar))
                    break
                except Exception:
                    time.sleep(0.02)
            if struct is None:
                raise RuntimeError("failed to read POSCAR")

            # run coarse relaxation
            t0 = time.time()
            calc.fmax = args.fmax_screen
            calc.max_steps = args.max_steps_screen
            res = calc.relax(struct)
            E_screen = float(res["energy"])
            struct_screen = res["final_structure"]
            elapsed = round(time.time() - t0, 3)

            # capacity check for waiting_pool
            ok = maybe_evict_for_capacity(POOL, args.pool_cap, E_screen)
            if not ok:
                # drop this candidate silently; pool is full with better structures
                sys.stderr.write(
                    f"[fast {k}] drop candidate E={E_screen:.6f} due to pool full.\n"
                )
            else:
                # insert into waiting_pool/<task_id> atomically
                task_id = f"{int(time.time() * 1e9)}_{k}_{uuid.uuid4().hex[:8]}"
                tmpd  = POOL / (".tmp_" + task_id)
                final = POOL / task_id
                tmpd.mkdir(parents=True, exist_ok=True)
                atomic_write_poscar(struct_screen, tmpd / "POSCAR")
                if savef.exists():
                    atomic_copy(savef, tmpd / "SAVE")
                meta = {
                    "task_id": task_id,
                    "source_slot": k,
                    "energy_screen": E_screen,
                    "fmax_screen": args.fmax_screen,
                    "max_steps_screen": args.max_steps_screen,
                    "elapsed_screen_s": elapsed,
                    "stamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                }
                atomic_write_json(tmpd / "meta.json", meta)
                os.rename(tmpd, final)

        except Exception as e:
            sys.stderr.write(f"[fast {k}] ERROR: {e}\n")

        finally:
            # tell C++ "this slot is free again"
            donef.touch()
            try:
                gof.unlink()
            except FileNotFoundError:
                pass

    sys.stderr.write(f"[fast {k}] bye\n")

if __name__ == "__main__":
    main()

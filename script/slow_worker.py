#!/usr/bin/env python3
"""
Slow worker: refine-only, consuming jobs from waiting_pool.
- Input pool:
    waiting_pool/<task_id>/{POSCAR, SAVE, meta.json}
- Each slow worker:
    * Scans waiting_pool, sorts candidates by 'energy_screen' ascending,
    * Atomically claims a job via rename:
          waiting_pool/<task_id> -> waiting_work/<task_id>
    * Runs a stricter relaxation (refinement),
    * Writes result into:
          refine_outbox/<task_id>/{CONTCAR, SAVE, energy, meta.json}
      where meta.json is updated with energy_final, timings, worker_id, etc.
    * Drops a lightweight report for the master into:
          reports/<task_id>.json
    * Removes the directory in waiting_work/<task_id> after finishing.
- Counters:
    ROOT/counters/slow_count incremented for each claimed job.
"""

import os, sys, json, time, tempfile, signal
from pathlib import Path
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

# ---------- waiting_pool selection ----------
def list_pool_sorted_by_energy(pool: Path):
    """Return list of task directories sorted by energy_screen ascending."""
    cands = []
    if pool.exists():
        for d in pool.iterdir():
            if not d.is_dir() or d.name.startswith(".tmp_"):
                continue
            try:
                meta = json.loads((d / "meta.json").read_text())
                E = float(meta.get("energy_screen", 1e99))
            except Exception:
                E = 1e99
            cands.append((E, d))
    cands.sort(key=lambda x: x[0])
    return [d for _, d in cands]

# ---------- main ----------
def main():
    import argparse
    ap = argparse.ArgumentParser(description="Slow worker (refine-only from waiting_pool)")
    ap.add_argument("--worker-id", type=str, default="slow-00")
    ap.add_argument("--root", type=str, default=".")
    ap.add_argument("--model", type=str, default="GRACE-2L-OMAT")
    ap.add_argument("--device", type=str, default="cpu", choices=["cpu", "cuda"])
    ap.add_argument("--fmax_refine", type=float, default=0.01)
    ap.add_argument("--max_steps_refine", type=int, default=400)
    ap.add_argument("--sleep_idle", type=float, default=0.2)
    args = ap.parse_args()

    ROOT = Path(args.root).resolve()
    POOL = ROOT / "waiting_pool"
    WORK = ROOT / "waiting_work"
    OUT  = ROOT / "refine_outbox"
    REPT = ROOT / "reports"

    for d in (WORK, OUT, REPT):
        d.mkdir(parents=True, exist_ok=True)

    calc = GraceCalculator(
        model=args.model,
        device=args.device,
        verbose=True,
        fmax=args.fmax_refine,
        max_steps=args.max_steps_refine,
    )
    sys.stderr.write(f"[{args.worker_id}] initialized @ {ROOT}\n")

    stop = {"flag": False}
    def _sigterm(*_): stop["flag"] = True
    signal.signal(signal.SIGINT, _sigterm)
    signal.signal(signal.SIGTERM, _sigterm)

    while not stop["flag"]:
        # Pick lowest-energy task from waiting_pool
        picked = None
        for d in list_pool_sorted_by_energy(POOL):
            target = WORK / d.name
            try:
                os.rename(d, target)  # atomic claim
                picked = target
                break
            except OSError:
                continue
        if picked is None:
            time.sleep(args.sleep_idle)
            continue

        try:
            atomic_increment_counter(ROOT, "slow_count")

            meta_in = {}
            meta_file = picked / "meta.json"
            try:
                meta_in = json.loads(meta_file.read_text())
            except Exception:
                meta_in = {}

            struct = Structure.from_file(str(picked / "POSCAR"))

            t0 = time.time()
            calc.fmax = args.fmax_refine
            calc.max_steps = args.max_steps_refine
            res = calc.relax(struct)
            E_final = float(res["energy"])
            struct_final = res["final_structure"]
            elapsed = round(time.time() - t0, 3)

            task_id = meta_in.get("task_id", picked.name)
            tmpd  = OUT / (".tmp_" + task_id)
            final = OUT / task_id
            tmpd.mkdir(parents=True, exist_ok=True)

            atomic_write_poscar(struct_final, tmpd / "CONTCAR")
            atomic_write_text(tmpd / "energy", f"{E_final:.12f}\n")
            meta_out = {
                **meta_in,
                "task_id": task_id,
                "energy_final": E_final,
                "fmax_refine": args.fmax_refine,
                "max_steps_refine": args.max_steps_refine,
                "elapsed_refine_s": elapsed,
                "worker_id": args.worker_id,
                "stamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            }
            atomic_write_json(tmpd / "meta.json", meta_out)

            # pass-through SAVE if present
            save_in = picked / "SAVE"
            if save_in.exists():
                atomic_copy(save_in, tmpd / "SAVE")

            os.rename(tmpd, final)

            # lightweight report for C++ master
            report = {
                "task_id": task_id,
                "energy_final": E_final,
                "energy_screen": meta_in.get("energy_screen"),
                "worker_id": args.worker_id,
                "stamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            }
            atomic_write_json(REPT / f"{task_id}.json", report)

        except Exception as e:
            task_id = meta_in.get("task_id", picked.name) if "meta_in" in locals() else picked.name
            atomic_write_json(REPT / f"{task_id}.json", {
                "task_id": task_id,
                "status": "error",
                "error": str(e),
                "worker_id": args.worker_id,
                "stamp": time.strftime("%Y-%m-%d %H:%M:%S"),
            })
        finally:
            # cleanup working dir
            try:
                for p in picked.iterdir():
                    try:
                        p.unlink()
                    except IsADirectoryError:
                        pass
                picked.rmdir()
            except Exception:
                pass

    sys.stderr.write(f"[{args.worker_id}] bye\n")

if __name__ == "__main__":
    main()

"""Rerun exact SCIP on every stored solver input listed in checks.tsv.

Usage (from the repository root):
    SCIP=/path/to/scip python3 code/rerun_exact.py logs/rerun --workers 16
Each case is solved with exact solving enabled before the problem is read. Only the
status 'problem is solved [infeasible]' together with 'exact/enable = TRUE' in the log
counts as an exclusion. A summary is written to <logdir>/summary.tsv.
"""
import argparse
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

root = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logdir", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--seconds", type=int, default=7200)
    parser.add_argument("--scip", default=os.environ.get("SCIP", "scip"))
    args = parser.parse_args()
    args.logdir.mkdir(parents=True, exist_ok=True)
    cases = [line.split("\t")[0] for line in (root / "checks.tsv").read_text().splitlines()[1:]]

    def run(case):
        model = (root / "models" / f"{case}.mps").resolve()
        command = [args.scip, "-c", "set exact enable TRUE",
                   "-c", f"set limits time {args.seconds}", "-c", f"read {model}",
                   "-c", "optimize", "-c", "quit"]
        job = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             stdin=subprocess.DEVNULL, text=True)
        (args.logdir / f"{case}.log").write_text(job.stdout)
        status = re.findall(r"SCIP Status\s*: (.*)", job.stdout)
        nodes = re.findall(r"Solving Nodes\s*: (\d+)", job.stdout)
        seconds = re.findall(r"Solving Time \(sec\)\s*: ([\d.]+)", job.stdout)
        exact = "exact/enable = TRUE" in job.stdout
        ok = exact and status == ["problem is solved [infeasible]"]
        row = [case, status[0] if status else "NO STATUS", nodes[0] if nodes else "",
               seconds[0] if seconds else "", "exact" if exact else "NOT EXACT",
               "OK" if ok else "UNRESOLVED"]
        print("\t".join(row), flush=True)
        return row

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(run, cases))
    header = ["case", "status", "nodes", "seconds", "mode", "verdict"]
    (args.logdir / "summary.tsv").write_text(
        "\n".join("\t".join(r) for r in [header] + rows) + "\n")
    bad = [r[0] for r in rows if r[-1] != "OK"]
    print("ALL CASES EXACT INFEASIBLE" if not bad else f"UNRESOLVED: {' '.join(bad)}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())

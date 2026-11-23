"""Run exact SCIP on individual cusp-order and diamond-gap cases.

Examples:
  python check_models.py 3
  python check_models.py 2 --gaps 1 2 3 --workers 3
  python check_models.py 0 --gaps 1 2 3 4 5 6 7 --nonzero-gap

The last command can require substantial time. Only an exact 'infeasible'
status is counted as an exclusion. A time limit is reported as unresolved.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import redirect_stdout
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory

from unit_model import write_model
from audit_mps import audit_mps


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("order", choices=range(5), type=int)
    parser.add_argument("--gaps", nargs="+", type=int)
    parser.add_argument("--nonzero-gap", action="store_true")
    parser.add_argument("--endpoint", type=int, choices=(-1,1))
    parser.add_argument("--prefix", type=int, nargs="+",
                        help="split by the first nonzero exponent in these positions")
    parser.add_argument("--zeros", type=int, nargs="+", default=[],
                        help="continue a previous prefix split in its all-zero branch")
    parser.add_argument("--reduced", action="store_true")
    parser.add_argument("--seconds", type=int, default=300)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--scip", default="scip")
    parser.add_argument("--logs", type=Path)
    args = parser.parse_args()
    if args.nonzero_gap and (args.order != 0 or not args.gaps):
        parser.error("--nonzero-gap requires order 0 and a gap")
    if args.endpoint is not None and not args.nonzero_gap:
        parser.error("--endpoint requires --nonzero-gap")
    gaps = args.gaps or [None]
    signs = [-1, 1] if args.nonzero_gap else [None]
    if args.endpoint is not None:
        signs = [args.endpoint]
    if (any(not 2 <= i <= 95 for i in args.zeros) or len(set(args.zeros)) != len(args.zeros)
            or set(args.zeros).intersection(args.prefix or [])):
        parser.error("zero positions must be valid, distinct and outside the new prefix")
    inherited = [(index, 0) for index in args.zeros]
    splits = [("", inherited)]
    if args.prefix:
        if len(set(args.prefix)) != len(args.prefix) or any(not 2 <= i <= 95 for i in args.prefix):
            parser.error("prefix positions must be distinct indices between 2 and 95")
        splits = []
        zeros = inherited.copy()
        for index in args.prefix:
            for value in (-1, 1):
                splits.append((f"_x{index}v{value}", zeros + [(index, value)]))
            zeros.append((index, 0))
        splits.insert(0, ("_zeros", zeros))
    tasks = []
    with TemporaryDirectory() as workspace:
        for gap in gaps:
            for sign in signs:
                for suffix, fixed in splits:
                    label = f"s{args.order}_g{gap}_e{sign}{suffix}"
                    directory = Path(workspace, label)
                    directory.mkdir()
                    with (directory / "model.lp").open("w") as f, redirect_stdout(f):
                        write_model(args.order, gap, sign, args.reduced, fixed)
                    subprocess.run(["glpsol", "--lp", "model.lp", "--wmps", "model.mps", "--check"],
                                   cwd=directory, stdout=subprocess.DEVNULL, check=True)
                    audit_mps(directory / "model.mps", args.order, gap, sign,
                              args.reduced, fixed)
                    tasks.append((label, directory))
        if args.logs:
            args.logs.mkdir(parents=True, exist_ok=True)

        def check(task):
            label, directory = task
            command = [args.scip, "-c", "set exact enable TRUE",
                       "-c", f"set limits time {args.seconds}", "-c", "read model.mps",
                       "-c", "optimize", "-c", "quit"]
            job = subprocess.run(command, cwd=directory, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, text=True, check=True)
            if args.logs:
                (args.logs / (label + ".log")).write_text(job.stdout)
            statuses = [line for line in job.stdout.splitlines() if "SCIP Status" in line]
            nodes = [line for line in job.stdout.splitlines() if "Solving Nodes" in line]
            if len(statuses) != 1 or "exact/enable = TRUE" not in job.stdout:
                raise RuntimeError(f"{label}: failed to read an exact solve status")
            print(label, statuses[0], *nodes, flush=True)
            return "problem is solved [infeasible]" in statuses[0]

        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            results = list(pool.map(check, tasks))
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())

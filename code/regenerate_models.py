"""Regenerate every solver input listed in checks.tsv and compare it with models/.

Usage (from the repository root):  python3 code/regenerate_models.py
For each case, the LP file is written by unit_model.write_model, converted to fixed
MPS by glpsol, audited row by row by audit_mps.audit_mps, and its SHA-256 is compared
with checks.tsv and with the stored file in models/.
"""
import hashlib
import subprocess
import sys
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parent))
from unit_model import write_model  # noqa: E402
from audit_mps import audit_mps  # noqa: E402

root = Path(__file__).resolve().parent.parent


def canonical(text):
    """Rename the rows of a fixed MPS file to R1, R2, ... in order of appearance.

    The stored inputs of nine cases were written by an earlier version of unit_model.py
    which left the endpoint row unnamed (glpsol then calls it R0000004). Row names do
    not affect the problem, so the comparison is made up to them as well.
    """
    lines = text.splitlines()
    names, section = {}, None
    for line in lines:
        if line and line[0] not in " *":
            section = line.split()[0]
        elif section == "ROWS" and line.startswith(" "):
            names[line.split()[1]] = f"R{len(names) + 1}"
    out = []
    for line in lines:
        tokens = line.split()
        if line.startswith(" ") and tokens:
            tokens = [names.get(t, t) for t in tokens]
            line = " ".join(tokens)
        out.append(line)
    return "\n".join(out)


def parse_fixed(text):
    if not text:
        return ()
    return tuple((int(a), int(v)) for a, v in (item.split("=") for item in text.split(",")))


def main():
    rows = [line.split("\t") for line in (root / "checks.tsv").read_text().splitlines()[1:]]
    failures = 0
    with TemporaryDirectory() as workspace:
        for case, order, gap, endpoint, fixed, _nodes, digest, status in rows:
            s = int(order)
            gap = int(gap) if gap else None
            endpoint = int(endpoint) if endpoint else None
            fixed = parse_fixed(fixed)
            directory = Path(workspace, case)
            directory.mkdir()
            with (directory / "model.lp").open("w") as f, redirect_stdout(f):
                write_model(s, gap, endpoint, False, fixed)
            subprocess.run(["glpsol", "--lp", "model.lp", "--wmps", "model.mps", "--check"],
                           cwd=directory, stdout=subprocess.DEVNULL, check=True)
            mps = directory / "model.mps"
            audit_mps(mps, s, gap, endpoint, False, fixed)
            stored_path = root / "models" / f"{case}.mps"
            audit_mps(stored_path, s, gap, endpoint, False, fixed)
            fresh = hashlib.sha256(mps.read_bytes()).hexdigest()
            stored = hashlib.sha256(stored_path.read_bytes()).hexdigest()
            if stored != digest:
                verdict, ok = "STORED FILE DOES NOT MATCH checks.tsv", False
            elif fresh == stored:
                verdict, ok = "identical", True
            elif canonical(mps.read_text()) == canonical(stored_path.read_text()):
                verdict, ok = "identical up to row names", True
            else:
                verdict, ok = "DIFFERENT", False
            failures += not ok
            print(case, "regenerated and audited (both files):", verdict, flush=True)
    print("ALL MODELS REPRODUCED" if failures == 0 else f"{failures} FAILURES")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

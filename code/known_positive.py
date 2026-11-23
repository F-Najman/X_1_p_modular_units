"""Test the integer programs on a known positive: the unit H_4 H_7 / (H_1 H_8) of degree 478.

Usage (from the repository root):  SCIP=/path/to/scip python3 code/known_positive.py
After the diamond <1/4> and inversion, H_4 H_7/(H_1 H_8) has exponents e_1 = 1, e_2 = -1,
e_48 = -1, e_46 = 1 and lies in the branch s0_g1_e-1_x48v-1 of checks.tsv. We take that
branch's model from unit_model.write_model, raise the degree bound 380 to 478, and check
with exact SCIP that
  (a) the model with all exponents fixed to this vector is feasible, with optimal value 478;
  (b) the model with the exponents free is feasible, with a solution of degree <= 478.
"""
import os
import re
import subprocess
import sys
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parent))
from unit_model import write_model  # noqa: E402

SCIP = os.environ.get("SCIP", "scip")
prefix = [3, 95, 94, 61, 41, 5, 8, 64, 4, 63, 38]
fixed = [(a, 0) for a in prefix] + [(48, -1)]
unit = {1: 1, 2: -1, 48: -1, 46: 1}


def model_text(bound, fix_all):
    buffer = StringIO()
    with redirect_stdout(buffer):
        write_model(0, 1, -1, False, fixed)
    lines = []
    for line in buffer.getvalue().splitlines():
        if line.startswith(" mass:"):
            line = re.sub(r"<= 380$", f"<= {bound}", line)
        line = re.sub(r"^ 0 <= (z\d+) <= 380$", rf" 0 <= \1 <= {bound}", line)
        match = re.match(r"^ -1 <= (e\d+) <= 1$", line)
        if fix_all and match:
            a = int(match.group(1)[1:])
            line = f" {match.group(1)} = {unit.get(a, 0)}"
        lines.append(line)
    text = "\n".join(lines) + "\n"
    assert text.count(f"<= {bound}") == 96 and "380" not in text
    return text


def solve(text):
    with TemporaryDirectory() as workspace:
        Path(workspace, "model.lp").write_text(text)
        subprocess.run(["glpsol", "--lp", "model.lp", "--wmps", "model.mps", "--check"],
                       cwd=workspace, stdout=subprocess.DEVNULL, check=True)
        command = [SCIP, "-c", "set exact enable TRUE", "-c", "set limits time 600",
                   "-c", "set limits solutions 1", "-c", "read model.mps",
                   "-c", "optimize", "-c", "display solution", "-c", "quit"]
        out = subprocess.run(command, cwd=workspace, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, text=True).stdout
    assert "exact/enable = TRUE" in out
    status = re.search(r"SCIP Status\s*: (.*)", out).group(1)
    value = re.search(r"objective value:\s*(\S+)", out)
    exps = {int(a): int(round(float(x))) for a, x in re.findall(r"^e(\d+)\s+(\S+)", out, re.M)}
    return status, value.group(1) if value else None, exps


status, value, exps = solve(model_text(478, True))
print("(a) fixed to H4H7/(H1H8):", status, "objective", value)
assert value is not None and int(value) == 478
status, value, exps = solve(model_text(478, False))
print("(b) exponents free:", status, "objective", value,
      "nonzero exponents", {a: x for a, x in sorted(exps.items()) if x})
assert value is not None and int(value) <= 478
print("Known positive successfully verified")

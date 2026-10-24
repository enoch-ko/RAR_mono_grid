#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""Isolated useful-bulgeless replacements. --plan only until separately approved.
Reuse a private module instance of the frozen, tested worker without editing it.
No submission/cancellation/refill code; one galaxy per exclusive shard claim.
"""
import importlib.util
import json
import re
from pathlib import Path
ROOT = Path(__file__).resolve().parent
PREP = ROOT / "outputs/grid_recovery_20261019/bulgeless47_preparation_v1"
MANIFEST = PREP / "manifest.json"
OUTPUT = ROOT / "outputs/grid_recovery_20261019/log_linear_bulgeless47_v1"
ENGINE = ROOT / "run_loglinear_recovery_20261019.py"
spec = importlib.util.spec_from_file_location("_isolated_bulgeless47_engine", ENGINE)
engine = importlib.util.module_from_spec(spec)
spec.loader.exec_module(engine)
SETTINGS = dict(engine.SETTINGS)
sha, write_new = engine.sha, engine.write_new
JOBS = {935526:5, 935528:7, 935529:8, 935531:10}

def validate(m):
    if m.get("schema") != "loglinear-bulgeless47-v1" or m.get("settings") != SETTINGS:
        raise ValueError("unexpected replacement schema/settings")
    if Path(m["output_root"]).resolve() != OUTPUT.resolve():
        raise ValueError("isolated output namespace required")
    for kind in ("source_sha256", "input_sha256"):
        if not m.get(kind):
            raise ValueError("missing dependency pins")
        for path, digest in m[kind].items():
            if not Path(path).is_file() or sha(path) != digest:
                raise ValueError("frozen dependency changed: " + path)
    audit_path = Path(m["allocation_audit"])
    if sha(audit_path) != m["allocation_audit_sha256"]:
        raise ValueError("allocation audit changed")
    audit = json.loads(audit_path.read_text())
    rows = audit["jobs"]
    if {int(j["job_id"]):int(j["shard"]) for j in rows} != JOBS:
        raise ValueError("original four-job allocation required")
    original = [g for j in rows for g in j["unique_bulgeless"]]
    if len(original) != 47 or len(set(original)) != 47:
        raise ValueError("47 unique original useful galaxies required")
    excluded = m.get("verified_completed_galaxies", [])
    # A newly completed old shard requires a fresh evidence review, not silent omission.
    if excluded:
        raise ValueError("completed-product exclusion requires a new reviewed manifest version")
    groups = m["shards"]
    names = [g for group in groups.values() for g in group]
    if (set(groups) != {f"{i:02d}" for i in range(1,48)} or
            any(len(group) != 1 for group in groups.values()) or
            len(names) != 47 or len(set(names)) != 47 or set(names) != set(original)):
        raise ValueError("exact 47 one-galaxy allocation required")
    if any(not re.fullmatch(r"[A-Za-z0-9_-]+", g) for g in names):
        raise ValueError("invalid galaxy name")
    norm = json.loads(Path(m["normalization"]).read_text())["galaxies"]
    if any(norm[g]["Lbulge_1e9"] != 0 or norm[g]["kappa"] != 1 for g in names):
        raise ValueError("all replacements must be zero-bulge kappa=1 controls")
    for key in (str(ENGINE), str(Path(__file__).resolve())):
        if key not in m["source_sha256"]:
            raise ValueError("launcher/engine pin missing")
    if str(Path(m["normalization"])) not in m["input_sha256"]:
        raise ValueError("normalization pin missing")
    return OUTPUT

def load_context():
    pin = json.loads((PREP / "launch_pin.json").read_text())
    if (pin["manifest_sha256"] != sha(MANIFEST) or
            pin["worker_sha256"] != sha(__file__) or pin["engine_sha256"] != sha(ENGINE)):
        raise ValueError("reviewed launch pin changed")
    m = json.loads(MANIFEST.read_text())
    validate(m)
    return m

# These assignments affect this private imported module only, never a running worker
# or the source/global state of the imported public 59-galaxy module in other processes.
engine.PREP, engine.MANIFEST = PREP, MANIFEST
engine.load_context = load_context
engine.validate = validate
build_plan, run = engine.build_plan, engine.run

def main():
    engine.main()

if __name__ == "__main__":
    main()

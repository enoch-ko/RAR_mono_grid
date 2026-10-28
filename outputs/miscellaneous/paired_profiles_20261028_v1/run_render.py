#!/usr/bin/env python3
"""Run the paired profile plotter with the current reorganized GP-results path."""
from pathlib import Path
import hashlib, json, sys, types

ROOT = Path.cwd()
sys.path.insert(0, str(ROOT))
PACKAGE = ROOT / "outputs/miscellaneous/paired_profiles_20261028_v1"
DRIVER = ROOT / "miscellaneous/support_code/plot_paired_forces_20261021.py"
PREP = PACKAGE / "manifest.json"
OUT = PACKAGE / "force_plots_v3"
EXPECTED_DRIVER_SHA = "ef06d1773bb4b45322497dedbbc1f662bea2117203320d8c49d230c2b3cbc5b6"
OLD_GP_PATH = "outputs/science_20261016/gp_d5000_a16"
CURRENT_GP_PATH = "outputs/GP_searches/science_20261016/gp_d5000_a16"

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

if sha(DRIVER) != EXPECTED_DRIVER_SHA:
    raise SystemExit("Canonical plot renderer hash changed; refusing to run")
source = DRIVER.read_text()
if source.count(OLD_GP_PATH) != 1:
    raise SystemExit("Expected exactly one legacy GP path in renderer")
source = source.replace(OLD_GP_PATH, CURRENT_GP_PATH, 1)
for galaxy in ("UGC04325", "KK98-251"):
    if not (ROOT / CURRENT_GP_PATH / f"{galaxy}_gp_search.json").is_file():
        raise SystemExit(f"Missing verified GP result for {galaxy}")

module = types.ModuleType("plot_paired_forces_20261021")
module.__file__ = str(DRIVER.resolve())
exec(compile(source, str(DRIVER), "exec"), module.__dict__)
module.ROOT = ROOT
module.PREP = PREP
module.OUT = OUT
args = sys.argv[1:]
if len(args) != 1 or args[0] not in ("--plan", "--run-forces"):
    raise SystemExit("Usage: run_render.py --plan | --run-forces")
sys.argv = [str(DRIVER), args[0]]
module.main()
if args[0] == "--run-forces":
    files = sorted(p for p in OUT.iterdir() if p.is_file())
    receipt = {
        "mode": args[0],
        "canonical_renderer": str(DRIVER),
        "canonical_renderer_sha256": sha(DRIVER),
        "runtime_renderer_sha256": hashlib.sha256(source.encode()).hexdigest(),
        "runtime_override": {"old_gp_path": OLD_GP_PATH, "current_gp_path": CURRENT_GP_PATH, "replacements": 1},
        "input_manifest_sha256": sha(PREP),
        "output_files_at_completion": {p.name: sha(p) for p in files},
    }
    (OUT / "runtime_receipt.json").write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")

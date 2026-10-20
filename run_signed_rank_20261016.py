#!/usr/bin/env python
"""Versioned signed-rank adapter around the unchanged normalization launcher.

PREPARED ONLY. Requires a separately reviewed/frozen signed-policy manifest and
NEW outputs/signed_rank_20261016 namespace. Never accepts the active manifest.
Use --describe for a no-science description. Queue approval remains separate.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
NEW_ROOT = ROOT / "outputs/signed_rank_20261016"
POLICY_ID = "finite_signed_gbar_rank_v1"
NEW_SOURCES = ("signed_rank_policy_20261016.py", "run_signed_rank_20261016.py")


def validate_manifest(manifest):
    if manifest.get("signed_rank_policy") != POLICY_ID:
        raise ValueError("new signed-rank manifest required; active manifests forbidden")
    if Path(manifest["pilot_output_root"]).resolve() != NEW_ROOT:
        raise ValueError("new, isolated signed-rank output namespace required")
    pins = manifest.get("source_sha256", {})
    for name in NEW_SOURCES:
        expected = pins.get(str(ROOT / name), pins.get(name))
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if expected != actual:
            raise ValueError("missing or mismatched new-source pin: " + name)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv == ["--describe"]:
        print(json.dumps({"policy": POLICY_ID, "status": "prepared_not_submitted",
            "required_new_output_root": str(NEW_ROOT),
            "required_new_sources": list(NEW_SOURCES),
            "science_evaluations": 0, "queue_submissions": 0}))
        return
    if not argv or argv[0] not in ("gp", "reference"):
        raise ValueError("mode must be gp/reference; no QP mode in this revision")
    if argv.count("--manifest") != 1 or "--" not in argv:
        raise ValueError("explicit manifest and forwarded argument separator required")
    prefix = argv[:argv.index("--")]
    manifest = json.loads(Path(prefix[prefix.index("--manifest") + 1]).read_text())
    validate_manifest(manifest)
    if any(x.startswith("--reference-json") or x.startswith("--gbar-floor-relative") for x in argv):
        raise ValueError("legacy QP import and applied gbar-floor options not supported")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    import signed_rank_policy_20261016 as policy
    policy.verify_parents()
    if argv[0] == "gp":
        policy.install_gp()
    else:
        policy.install_reference()
    import run_corrected_recovery_pilot as recovery
    original = sys.argv
    try:
        sys.argv = [str(ROOT / "run_corrected_recovery_pilot.py")] + argv
        recovery.main()
    finally:
        sys.argv = original


if __name__ == "__main__":
    main()

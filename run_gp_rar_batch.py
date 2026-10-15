#!/usr/bin/env python
"""Prepare or execute 25 frozen seven-galaxy GP shards; never submit jobs."""
import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"
import argparse
import collections
import hashlib
import json
import shlex
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT/"outputs/gp_rar_full_20261013_d5000_a16"
SETTINGS = {"directions": 5000, "amplitude_steps": 16, "ell_over_rd": [.5, 1., 2., 4.],
            "seed": 20261013, "replicate": 0, "sampling_scheme": "nested-v1",
            "checkpoint_mode": "compact", "min_ml": .01, "max_ml": 10.,
            "length_fraction": .3, "gbar_floor_relative": 1e-8,
            "gas_ml": 1., "bulge_ml": .7, "gas_ml_dex": .04, "bulge_ml_dex": .1,
            "ml_samples": 9, "prob_grid_size": 400, "zpadding": 8.,
            "boundary_log_tolerance": 2e-5, "boundary_max_refinements": 8}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def allocate(items, nshards=25, capacity=7):
    if len(items) != nshards*capacity or len({x["galaxy"] for x in items}) != len(items):
        raise ValueError("allocation requires exactly 175 unique galaxies")
    shards = [{"shard": i+1, "galaxies": [], "estimated_weight": 0} for i in range(nshards)]
    for item in sorted(items, key=lambda x: (-x["weight"], x["galaxy"])):
        shard = min((s for s in shards if len(s["galaxies"]) < capacity),
                    key=lambda s: (s["estimated_weight"], s["shard"]))
        shard["galaxies"].append(item["galaxy"])
        shard["estimated_weight"] += item["weight"]
    for shard in shards:
        shard["galaxies"].sort()
    return shards


def payload(galaxy, output, reference=None):
    command = [str(ROOT/"gp_rar_search.py"), "--galaxy", galaxy,
               "--output-root", str(output/"galaxies")]
    for key, value in SETTINGS.items():
        command.append("--"+key.replace("_", "-"))
        command.extend(str(v) for v in value) if isinstance(value, list) else command.append(str(value))
    if reference:
        command.extend(["--reference-json", galaxy+"="+reference["path"]])
    return command


def queue_commands(output):
    base = 'addqueue -q cmb -n 1 -m 8 -c "5 hrs" '
    return [base+shlex.join([str(ROOT/"run_gp_rar_batch.py"), "--output-root", str(output),
                             "--shard", str(i)]) for i in range(1, 26)]


def prepare(output):
    import numpy as np
    import gp_rar_search as gp
    import inverse_rar_diagnostic as inverse
    from utils import load_sparc_table, DEFAULT_DENS_ROOT, DEFAULT_TABLE_PATH
    if output.exists():
        raise FileExistsError("choose a new output root")
    table = load_sparc_table()
    names = table["Galaxy"].tolist()
    if len(names) != 175 or len(set(names)) != 175:
        raise ValueError("expected 175 unique SPARC galaxies")
    sources = {str(ROOT/name): sha(ROOT/name) for name in
               (*gp.SOURCE_FILES, "gp_rar_search.py", "run_gp_rar_batch.py")}
    validation = sorted((ROOT/"outputs").glob(inverse.VALIDATION_GLOB))
    sources.update({str(p): sha(p) for p in validation})
    input_paths = {str(Path(DEFAULT_TABLE_PATH)): sha(DEFAULT_TABLE_PATH)}
    items, references, reference_notes = [], {}, {}
    args = SimpleNamespace(**SETTINGS)
    args.length_fraction = SETTINGS["length_fraction"]
    for galaxy in names:
        gate = inverse.cache_gate(galaxy, ROOT/"outputs")
        dens = Path(DEFAULT_DENS_ROOT)/f"{galaxy}.dens"
        rotmod = Path(gp.mono.DEFAULT_DATA_ROOT)/f"{galaxy}_rotmod.dat"
        hashes = {"dens": sha(dens), "rotmod": sha(rotmod), "table": sha(DEFAULT_TABLE_PATH)}
        input_paths.update({str(dens): hashes["dens"], str(rotmod): hashes["rotmod"]})
        rad, brightness, r, vobs, errv, vgas, vbul, zd, excluded = inverse.galaxy_data(
            galaxy, gp.mono.DEFAULT_DATA_ROOT, DEFAULT_DENS_ROOT)
        if len(rad) < 2 or len(r) < 2:
            raise ValueError(f"insufficient retained grid: {galaxy}")
        rd = float(table.loc[table.Galaxy == galaxy].iloc[0]["Rdisk"])
        if not np.isfinite(rd) or rd <= 0:
            raise ValueError(f"invalid Rdisk: {galaxy}")
        # Proxy balances knot matrix work AND probability track size, not a runtime guarantee.
        items.append({"galaxy": galaxy, "knots": len(rad), "retained_points": len(r),
                      "weight": int(len(rad)**2*len(r)), "cache_gate": gate})
        candidates = [ROOT/"outputs/rerun_boundary_20261006/binned_ML_opt/galaxies"/f"{galaxy}_inverse.json"]
        if galaxy == "UGC02916":
            candidates.insert(0, ROOT/"outputs/binned_ML_opt/osqp_UGC02916_tight_20261006/UGC02916_warm_start.json")
        reasons = []
        for path in candidates:
            if not path.is_file():
                reasons.append("missing "+str(path))
                continue
            try:
                xi, info = gp.reference_profile(path, galaxy, np.asarray(rad), .01, 10., args, hashes)
            except (ValueError, KeyError, TypeError) as error:
                reasons.append(str(error))
                continue
            references[galaxy] = {"path": str(path), "sha256": sha(path), "audit": info}
            break
        reference_notes[galaxy] = {"selected": galaxy in references, "rejected_reasons": reasons}
    shards = allocate(items)
    manifest = {"status": "prepared_not_submitted", "settings": SETTINGS,
                "output_root": str(output), "source_and_validation_sha256": sources,
                "input_sha256": input_paths, "galaxies": items, "shards": shards,
                "references": references, "reference_notes": reference_notes,
                "sampling": "per-galaxy/per-ell SHA256 streams, replicate0, prefix-nested budgets",
                "existing_pilots": "preserved separately; pool only after explicit compatibility audit",
                "compact_checkpoint_policy": "all final xi in indexed profiles.npy; partial progress per ell; no resume",
                "interpretation": "finite-search first survey, not converged/global optima"}
    for path, expected in sources.items():
        if sha(path) != expected:
            raise RuntimeError("source changed during preparation")
    output.mkdir(parents=True)
    write_json(output/"batch_manifest.json", manifest)
    review = ["# GP full-sample submission review", "", "PREPARED ONLY. 25 jobs x 7 galaxies = 175.",
              "1 CPU and 8 GB per job; up to 25 CPUs / 200 GB concurrently.",
              "Scheduling comment 5 hrs is not a wall-time limit. Runtime is galaxy-dependent.",
              "5000 directions per each of 4 ell/Rd, 16 amplitudes, unchanged physics/probability controls.",
              "Compact checkpoints; shared galaxies/ outputs, per-shard reports and per-galaxy logs.",
              "First finite-search survey; frontier convergence remains unestablished.", ""]
    review = review+[chr(96)*3+"bash"]+queue_commands(output)+[chr(96)*3, "", "## Galaxy allocation", ""]
    review += [f"Shard {s['shard']:02d}: "+", ".join(s["galaxies"]) for s in shards]
    review += ["", "Read-only preflight: python3 -B run_gp_rar_batch.py --dry-run", "",
               "After completion: python3 -B run_gp_rar_batch.py --audit", "No CSV merger is needed for these JSON/NPY products."]
    with (output/"submission_review.md").open("x") as stream:
        stream.write("\n".join(review)+"\n")
    print("PREPARED", output, "galaxies", len(items), "references", len(references), flush=True)
    print("SHARD_SIZES", [len(s["galaxies"]) for s in shards])
    print("NO_JOBS_SUBMITTED")


def verify_manifest(output):
    manifest = json.loads((output/"batch_manifest.json").read_text())
    if manifest["settings"] != SETTINGS or manifest["output_root"] != str(output):
        raise ValueError("batch settings/output mismatch")
    for path, expected in manifest["source_and_validation_sha256"].items():
        if sha(path) != expected:
            raise RuntimeError("source/validation changed: "+path)
    for path, expected in manifest["input_sha256"].items():
        if sha(path) != expected:
            raise RuntimeError("input changed: "+path)
    for ref in manifest["references"].values():
        if sha(ref["path"]) != ref["sha256"]:
            raise RuntimeError("reference changed: "+ref["path"])
    return manifest


def run_shard(output, shard_id):
    manifest = verify_manifest(output)
    shard = manifest["shards"][shard_id-1]
    report = output/f"shard_{shard_id:02d}.json"
    if report.exists():
        raise FileExistsError(report)
    (output/"logs").mkdir(exist_ok=True)
    results = []
    for galaxy in shard["galaxies"]:
        verify_manifest(output)
        command = payload(galaxy, output, manifest["references"].get(galaxy))
        print("START", galaxy, shlex.join(command), flush=True)
        started = time.monotonic()
        log = output/"logs"/f"{galaxy}.log"
        with log.open("x") as stream:
            done = subprocess.run([sys.executable, "-B", *command], cwd=ROOT,
                                  stdout=stream, stderr=subprocess.STDOUT)
        entry = {"galaxy": galaxy, "returncode": done.returncode,
                 "seconds": time.monotonic()-started, "log": str(log)}
        if done.returncode == 0:
            path = output/"galaxies"/f"{galaxy}_gp_search.json"
            with path.open() as stream:
                data = json.load(stream)
            numerical = sum(not r.get("valid") and r.get("reason") not in
                            ("signed_gbar_below_floor", "nonfinite_signed_gbar") for r in data["records"])
            nonfinite_fields = sum(not r.get("valid") and r.get("reason") == "nonfinite_signed_gbar"
                                   for r in data["records"])
            entry.update(status=data["status"], records=len(data["records"]),
                         numerical_failures=numerical+nonfinite_fields)
            del data
        else:
            entry["status"] = "process_failure"
        results.append(entry)
        write_json(output/f"{galaxy}_completion.json", entry)
        print("FINISH", galaxy, json.dumps(entry), flush=True)
    verify_manifest(output)
    failed = sum(x.get("status") != "finite_search_completed" or x.get("numerical_failures", 0)>0
                 for x in results)
    write_json(report, {"shard": shard_id, "expected": shard["galaxies"],
                        "results": results, "failed_or_numerically_invalid": failed,
                        "manifest_sha256": sha(output/"batch_manifest.json")})
    if failed:
        raise SystemExit(1)


def audit(output):
    manifest = verify_manifest(output)
    names = [x["galaxy"] for x in manifest["galaxies"]]
    missing, statuses = [], collections.Counter()
    for galaxy in names:
        path = output/f"{galaxy}_completion.json"
        if not path.exists():
            missing.append(galaxy)
            continue
        data = json.loads(path.read_text())
        statuses[data["status"]] += 1
    print(json.dumps({"expected": len(names), "completion_statuses": dict(statuses),
                      "missing": missing, "completed_shard_reports": len(list(output.glob("shard_*.json"))),
                      "scientific_validation": "separate per-result QC still required"}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--shard", type=int, choices=range(1, 26))
    mode.add_argument("--audit", action="store_true")
    args = parser.parse_args()
    output = args.output_root.resolve()
    if args.prepare:
        prepare(output)
    elif args.dry_run:
        manifest = verify_manifest(output)
        names = [g for s in manifest["shards"] for g in s["galaxies"]]
        assert len(names) == len(set(names)) == 175
        assert all(len(s["galaxies"]) == 7 for s in manifest["shards"])
        for galaxy in names:
            if (output/"galaxies"/f"{galaxy}_gp_search.json").exists() or (output/"galaxies"/f"{galaxy}_gp_checkpoints").exists():
                raise FileExistsError("existing galaxy output: "+galaxy)
        print(json.dumps({"settings": SETTINGS, "shards": manifest["shards"],
                          "commands": queue_commands(output), "jobs_submitted": 0}, indent=2))
    elif args.audit:
        audit(output)
    else:
        run_shard(output, args.shard)


if __name__ == "__main__":
    main()

#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""Isolated, fixed-allocation dispatcher for signed-rank GP/reference shards (2026-10-19).
No retries/refills; dry-run is the default. Manifest is immutable and versioned.
"""
import argparse, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
POLICY = "finite_signed_gbar_rank_v1"
SCHEMA = "signed-rank-shards-v1"
OUTPUT_PARENT = (ROOT / "outputs/signed_rank_20261016").resolve()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def inside(path, parent):
    try:
        Path(path).resolve().relative_to(Path(parent).resolve())
        return True
    except ValueError:
        return False


def validate_manifest(m, manifest_path):
    if m.get("schema") != SCHEMA or m.get("signed_rank_policy") != POLICY:
        raise ValueError("unsupported manifest schema or signed-rank policy")
    out = Path(m["output_root"]).resolve()
    runner_manifest = Path(m["runner_manifest"]).resolve()
    if not runner_manifest.is_file() or sha256(runner_manifest) != m.get("runner_manifest_sha256"):
        raise ValueError("runner manifest missing or hash mismatch")
    if not inside(out, OUTPUT_PARENT) or out == OUTPUT_PARENT:
        raise ValueError("output_root must be a new isolated child of signed_rank_20261016")
    if "integration_smoke_v1_20261016" in out.parts or "science_20261016" in out.parts:
        raise ValueError("pilot and production evidence roots are immutable")
    pins = m.get("source_sha256", {})
    if not pins:
        raise ValueError("manifest must pin dispatcher and all run sources")
    for name, expected in pins.items():
        p = Path(name)
        if not p.is_absolute():
            p = ROOT / p
        if not p.is_file() or sha256(p) != expected:
            raise ValueError("missing or changed source pin: " + str(name))
    shards = m.get("shards", {})
    if not shards or any(not isinstance(k, str) or not v for k, v in shards.items()):
        raise ValueError("manifest must contain non-empty fixed shard allocations")
    galaxies = [g for values in shards.values() for g in values]
    if len(galaxies) != len(set(galaxies)) or any(not isinstance(g, str) or not g for g in galaxies):
        raise ValueError("galaxy allocation contains duplicates or invalid names")
    mode = m.get("mode", "gp")
    if mode not in ("gp", "reference"):
        raise ValueError("unsupported component mode")
    if mode == "reference" and (set(m.get("seeds", {})) != set(galaxies) or
            any(type(s) is not int or s < 0 for s in m["seeds"].values())):
        raise ValueError("reference seeds must cover the exact allocation")
    template = m.get("argv_template")
    if not isinstance(template, list) or not template or any(not isinstance(x, str) for x in template):
        raise ValueError("argv_template must be a non-empty argument list, never a shell string")
    forbidden = ("--reference-json", "--gbar-floor-relative", "--replace-existing", "qp")
    if any(any(bad in token for bad in forbidden) for token in template):
        raise ValueError("legacy floor/reference import, QP, or overwrite option is forbidden")
    if m.get("max_runtime_seconds", 0) < 1 or m.get("minimum_free_bytes", 0) < 0:
        raise ValueError("runtime and storage guards are required")
    # Caller may only invoke this dispatcher with the exact same frozen manifest.
    if not Path(manifest_path).is_file():
        raise ValueError("manifest file is missing")
    return out


def render(template, galaxy, output_root, manifest_path, seed=None):
    fields = {"galaxy": galaxy, "output_root": str(output_root / galaxy),
              "manifest": str(Path(manifest_path).resolve()), "python": sys.executable,
              "seed": str(seed), "output_file": str(output_root / galaxy / (galaxy + "_B1024.json"))}
    return [token.format(**fields) for token in template]


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp-" + str(os.getpid()))
    tmp.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def dispatch(m, manifest_path, shard_id, dry_run=True):
    out = validate_manifest(m, manifest_path)
    if shard_id not in m["shards"]:
        raise ValueError("unknown shard id; dispatcher never refills another shard")
    if dry_run:
        return [render(m["argv_template"], g, out, m["runner_manifest"], m.get("seeds", {}).get(g)) for g in m["shards"][shard_id]]
    out.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(out).free < m["minimum_free_bytes"]:
        raise OSError("storage guard refused shard before launching any galaxy")
    commands = []
    for galaxy in m["shards"][shard_id]:
        target = out / galaxy
        if target.exists():
            raise FileExistsError("refusing existing per-galaxy output directory: " + galaxy)
        suffix = "_B1024.json" if m.get("mode") == "reference" else "_gp_search.json"
        result = target / (galaxy + suffix)
        status = target / "dispatcher_status.json"
        claim = out / ".claims" / (shard_id + "__" + galaxy + ".lock")
        claim.parent.mkdir(parents=True, exist_ok=True)
        try:
            fd = os.open(claim, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        except FileExistsError as exc:
            raise RuntimeError("already claimed; no automatic retry: " + galaxy) from exc
        os.close(fd)
        target.mkdir(parents=True, exist_ok=False)
        if result.exists() or status.exists():
            raise FileExistsError("refusing to overwrite existing per-galaxy output/status: " + galaxy)
        argv = render(m["argv_template"], galaxy, out, m["runner_manifest"], m.get("seeds", {}).get(galaxy))
        record = {"galaxy": galaxy, "shard": shard_id, "argv": argv,
                  "started_unix": time.time(), "state": "running"}
        write_json(status, record)
        env = os.environ.copy()
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            env[key] = "1"
        try:
            proc = subprocess.run(argv, cwd=ROOT, env=env, capture_output=True, text=True,
                                  timeout=m["max_runtime_seconds"], check=False)
            output_missing = proc.returncode == 0 and not result.is_file()
            record.update(state="complete" if proc.returncode == 0 and not output_missing else "failed",
                          returncode=proc.returncode, stdout=proc.stdout[-12000:],
                          stderr=proc.stderr[-12000:], finished_unix=time.time())
            if output_missing:
                record["error"] = "runner returned success but expected per-galaxy JSON is missing: " + str(result)
        except Exception as exc:
            record.update(state="failed", error=repr(exc), finished_unix=time.time())
            write_json(status, record)
            raise
        write_json(status, record)
        commands.append(record)
        if record["state"] != "complete":
            raise RuntimeError("galaxy failed; stopping shard without retry/refill: " + galaxy)
    return commands


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("batch", nargs="?", choices=("gp", "reference"))
    p.add_argument("shard_id", nargs="?")
    p.add_argument("--manifest", type=Path)
    p.add_argument("--shard")
    p.add_argument("--execute", action="store_true", help="otherwise print commands only")
    a = p.parse_args()
    if a.batch:
        if a.manifest or a.shard or not a.shard_id:
            p.error("positional batch requires a shard ID and no manifest override")
        if not os.environ.get("SLURM_JOB_ID"):
            p.error("queued positional execution requires SLURM_JOB_ID")
        a.manifest = ROOT / "outputs/signed_rank_20261016/preparation/full_20261019_v1" / (a.batch + "_dispatch.json")
        a.shard, a.execute = a.shard_id, True
    if not a.manifest or not a.shard:
        p.error("explicit manifest/shard required for dry-run")
    m = json.loads(a.manifest.read_text())
    result = dispatch(m, a.manifest, a.shard, dry_run=not a.execute)
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()


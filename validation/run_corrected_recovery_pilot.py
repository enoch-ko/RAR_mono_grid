#!/mnt/users/koe/RAR_mono_grid/.venv-osqp/bin/python
"""Opt-in corrected+normalized pilot launcher. No queue actions or legacy resume.
A frozen manifest is mandatory. Scientific modes require separate run approval.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import sys

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(manifest):
    for field in ("source_sha256", "input_sha256"):
        for path, expected in manifest[field].items():
            if sha(path) != expected:
                raise ValueError("frozen dependency changed: " + path)


def normalization_proof(manifest, manifest_path):
    return dict(policy="disc_to_catalogue_total_minus_published_bulge_v1",
                normalization_sha256=sha(manifest["normalization"]),
                recovery_manifest_sha256=sha(manifest_path),
                recovery_launcher_sha256=sha(__file__))


def install(manifest, proof):
    """Process-local correction only; never edit frozen source modules."""
    import numpy as np
    import utils
    import mono_SPARC_RAR as mono
    import inverse_rar_diagnostic_corrected as inverse
    norm = json.loads(Path(manifest["normalization"]).read_text())
    if norm.get("policy") != proof["policy"]:
        raise ValueError("unsupported normalization policy")
    records = norm["galaxies"]
    original = utils.get_SBdisk

    def normalized(galaxy, dens_root=utils.DEFAULT_DENS_ROOT, table=None):
        rad, sb = original(galaxy, dens_root=dens_root, table=table)
        rec = records[galaxy]
        cat = (utils.load_sparc_table() if table is None else table)
        row = cat.loc[cat.Galaxy == galaxy]
        if len(row) != 1 or not np.isclose(float(row.iloc[0].L), rec["Ltotal_1e9"], rtol=0, atol=1e-10):
            raise ValueError("normalization catalogue mismatch")
        if not np.isclose(rec["Ldisk_1e9"] + rec["Lbulge_1e9"], rec["Ltotal_1e9"], rtol=1e-12):
            raise ValueError("component luminosities do not sum to total")
        if not np.isclose(rec["kappa"], rec["Ldisk_1e9"] / rec["Ltotal_1e9"], rtol=1e-12):
            raise ValueError("invalid disc luminosity scale")
        corrected = sb * rec["kappa"]
        if not np.isclose(utils.compute_A_B(rad, corrected)[0], rec["Ldisk_1e9"] * 1e9, rtol=1e-12):
            raise ValueError("normalized luminosity integral mismatch")
        return rad, corrected

    utils.get_SBdisk = normalized
    mono.get_SBdisk = normalized
    inverse.get_SBdisk = normalized
    old_gate = inverse.cache_gate

    def gate(galaxy, root):
        evidence = old_gate(galaxy, root)
        return dict(evidence, normalization=proof,
                    normalization_qc="historical basis gate; current normalized direct/cache checks required")

    inverse.cache_gate = gate
    return inverse, records


def main():
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("mode", choices=("qp", "gp", "reference"))
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--seconds", type=int, default=900)
    # Module arguments are placed after an explicit -- to prevent ambiguity.
    prefix, forwarded = sys.argv[1:], []
    if "--" in prefix:
        index = prefix.index("--")
        prefix, forwarded = prefix[:index], prefix[index+1:]
    args = parser.parse_args(prefix)
    if args.seconds < 1 or not forwarded:
        parser.error("positive time budget and forwarded module arguments required")
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema") != "corrected-recovery-pilot-v1":
        raise ValueError("incompatible recovery manifest")
    verify(manifest)
    if "--replace-existing" in forwarded:
        raise ValueError("pilot outputs must not replace existing scientific products")
    output_option = "--output-root" if args.mode == "gp" else "--output"
    if output_option not in forwarded:
        raise ValueError("explicit new pilot output path required")
    target = Path(forwarded[forwarded.index(output_option)+1]).resolve()
    pilot_root = Path(manifest["pilot_output_root"]).resolve()
    if pilot_root not in target.parents:
        raise ValueError("output must be below frozen pilot namespace")
    for option, expected in (("--data-root", manifest["data_root"]),
                             ("--dens-root", manifest["dens_root"]),
                             ("--table", manifest["table"]),
                             ("--validation-root", manifest["validation_root"])):
        if any(v.startswith(option+"=") for v in forwarded):
            raise ValueError("use separate option/value tokens: " + option)
        if option in forwarded and Path(forwarded[forwarded.index(option)+1]).resolve() != Path(expected).resolve():
            raise ValueError("input root differs from frozen manifest: " + option)
    proof = normalization_proof(manifest, args.manifest)
    print("FROZEN_CORRECTED_RECOVERY_VERIFIED", json.dumps(proof), flush=True)
    # Bound native BLAS thread pools before importing any scientific modules.
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"

    def timed_out(*unused):
        raise TimeoutError("declared corrected-recovery pilot budget exhausted")

    signal.signal(signal.SIGALRM, timed_out)
    signal.alarm(args.seconds)
    try:
        if args.mode == "reference":
            # The reference driver applies kappa itself; do NOT install a second scaling.
            import common_reference_driver as driver
            if "--normalization" not in forwarded:
                forwarded += ["--normalization", manifest["normalization"]]
            elif Path(forwarded[forwarded.index("--normalization")+1]).resolve() != Path(manifest["normalization"]).resolve():
                raise ValueError("reference normalization differs from frozen configuration")
            sys.argv = [driver.__file__] + forwarded
            driver.main()
        else:
            inverse, records = install(manifest, proof)
            if args.mode == "qp":
                import inverse_rar_boundary_corrected as driver
                analyze = inverse.analyze

                def analyzed(options):
                    result = analyze(options)
                    result["normalization"] = dict(proof, galaxy=options.galaxy, **records[options.galaxy])
                    result["provenance"]["recovery_configuration"] = proof
                    return result

                inverse.analyze = analyzed
                sys.argv = [driver.__file__] + forwarded
                driver.main()
            else:
                import gp_rar_search_retained_corrected as driver
                for i, value in enumerate(forwarded):
                    if value == "--reference-json":
                        galaxy, marker, path = forwarded[i+1].partition("=")
                        saved = json.loads(Path(path).read_text())
                        if not marker or saved.get("provenance", {}).get("recovery_configuration") != proof:
                            raise ValueError("QP reference lacks this integrated normalization/configuration")
                        if saved.get("galaxy") != galaxy:
                            raise ValueError("reference galaxy mismatch")
                    if value.startswith("--reference-json="):
                        raise ValueError("use separate --reference-json GALAXY=PATH tokens")
                driver.SOURCE_FILES = tuple(dict.fromkeys((*driver.SOURCE_FILES,
                    Path(__file__).name, manifest["normalization"], str(args.manifest.resolve()))))
                run = driver.galaxy_run

                def galaxy_run(galaxy, options, reference=None):
                    result = run(galaxy, options, reference)
                    result["normalization"] = dict(proof, galaxy=galaxy, **records[galaxy])
                    result.setdefault("provenance", {})["recovery_configuration"] = proof
                    return result

                driver.galaxy_run = galaxy_run
                sys.argv = [driver.__file__] + forwarded
                driver.main()
        verify(manifest)
    finally:
        signal.alarm(0)


if __name__ == "__main__":
    main()

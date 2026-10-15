#!/usr/bin/env python
"""Isolated cached-field grid entry point with boundary-resolving probabilities.
Existing shared solvers and active legacy jobs are not changed.
"""
import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"
import argparse
import hashlib
import json
import sys
from pathlib import Path
from functools import partial
import mono_SPARC_RAR as mono
import probability_recursion as probability

HERE = Path(__file__).resolve().parent

def hashes():
    return {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
            for name in ("mono_SPARC_RAR.py", "utils.py", "probability_recursion.py",
                         "mono_SPARC_RAR_boundary.py")}

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--boundary-log-tolerance", type=float, default=2e-5)
    parser.add_argument("--boundary-max-refinements", type=int, default=8)
    parser.add_argument("--dry-run", action="store_true")
    options, rest = parser.parse_known_args()
    if not 0 < options.boundary_log_tolerance < float("inf") or options.boundary_max_refinements < 1:
        parser.error("finite positive tolerance and positive refinement count required")
    sys.argv = [str(HERE / "mono_SPARC_RAR.py")] + rest
    args = mono.parse_args()
    args.probability_method = "boundary_log_cumulative_v1"
    args.boundary_log_tolerance = options.boundary_log_tolerance
    args.boundary_max_refinements = options.boundary_max_refinements
    sources = hashes()
    args.numerical_source_sha256 = json.dumps(sources, sort_keys=True)
    if options.dry_run:
        print(json.dumps({"arguments": {k: str(v) for k, v in vars(args).items()},
                          "source_sha256": sources}, indent=2))
        return
    # Per-process injection only; the original nuisance integration, signed
    # field cache, physical prior mass accounting and CSV schemas are retained.
    mono.parse_args = lambda: args
    mono.log_monotonic_probability_for_gbar = partial(
        probability.boundary_log_probability_for_gbar,
        log_tolerance=options.boundary_log_tolerance,
        max_refinements=options.boundary_max_refinements)
    mono.source_sha256 = lambda: hashlib.sha256(
        json.dumps(sources, sort_keys=True).encode()).hexdigest()
    original = mono.process_galaxy
    def guarded(*positional, **keywords):
        if hashes() != sources:
            raise RuntimeError("Source changed during shard; refusing mixed algorithms")
        result = original(*positional, **keywords)
        if hashes() != sources:
            raise RuntimeError("Source changed within galaxy; refusing mixed provenance")
        return result
    mono.process_galaxy = guarded
    mono.main()

if __name__ == "__main__":
    main()


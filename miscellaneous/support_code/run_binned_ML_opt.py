#!/usr/bin/env python
"""Descriptive entry point for the positive-net-gbar binned M/L optimizer.
Use --dry-run for allocation and --audit after all ten shards finish.
"""
import os
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"
import sys
import run_inverse_rar_batch as batch

if __name__ == "__main__":
    batch.DEFAULT_ROOT = batch.HERE / "outputs/binned_ML_opt"
    # Include launcher source in shard guards, preserving the established
    # galaxy provenance and backward-compatible inverse result filenames.
    original_hashes = batch.source_hashes
    def hashes():
        return dict(original_hashes(), binned_ML_opt_launcher=batch.digest(__file__))
    batch.source_hashes = hashes
    sys.exit(batch.main())


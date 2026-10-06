#!/usr/bin/env python3
"""Validate and combine independently written SPARC run shards.

Each analysis job writes uniquely named shard CSVs and metadata into the shared
profile output directory. Run this merger after all jobs in that batch finish.
"""

import argparse
import csv
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path


FINAL_NAMES = {
    "profile_grid": "profile_grid.csv",
    "fully_marginalised": "fully_marginalised.csv",
    "metadata": "run_metadata.json",
}
SHARED_ARGUMENT_EXCLUSIONS = {
    "galaxy", "galaxy_range", "overwrite_shard", "shard_id", "testing",
}
PROFILE_GRID_REDUNDANT_COLUMNS = {
    "gas_ml_median", "disk_ml_median", "bulge_ml_median",
    "gas_ml_dex", "disk_ml_dex", "bulge_ml_dex",
    "ml_samples_per_term", "gas_term_convention", "tie_atol",
}
FULLY_MARGINALISED_REDUNDANT_COLUMNS = {
    "disk_ml_median", "disk_ml_dex", "slope_min", "slope_max", "grid_size",
}


def compact_fields(fieldnames, redundant_columns):
    """Return a CSV schema without run-level settings repeated on every row."""
    return [name for name in fieldnames if name not in redundant_columns]


def compact_rows(rows, redundant_columns):
    """Drop run-level columns from in-memory rows without editing source shards."""
    return [
        {name: value for name, value in row.items() if name not in redundant_columns}
        for row in rows
    ]


def read_csv(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            return None, []
        return list(reader.fieldnames), list(reader)


def read_shard(output_root, shard_id):
    suffix = f".shard_{shard_id:02d}"
    paths = {
        "profile_grid": output_root / f"profile_grid{suffix}.csv",
        "fully_marginalised": output_root / f"fully_marginalised{suffix}.csv",
        "metadata": output_root / f"run_metadata{suffix}.json",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Shard {shard_id} is incomplete; missing: {', '.join(missing)}")
    with paths["metadata"].open() as stream:
        metadata = json.load(stream)
    grid_fields, grid_rows = read_csv(paths["profile_grid"])
    full_fields, full_rows = read_csv(paths["fully_marginalised"])
    if grid_fields is None:
        raise ValueError(f"Shard {shard_id} has an empty profile-grid CSV")
    actual_shard_id = metadata.get("arguments", {}).get("shard_id")
    if str(actual_shard_id) != str(shard_id):
        raise ValueError(
            f"Shard file {shard_id} contains metadata for shard {actual_shard_id!r}"
        )
    return {
        "id": shard_id,
        "paths": paths,
        "metadata": metadata,
        "grid_fields": grid_fields,
        "grid_rows": compact_rows(grid_rows, PROFILE_GRID_REDUNDANT_COLUMNS),
        "full_fields": full_fields,
        "full_rows": compact_rows(full_rows, FULLY_MARGINALISED_REDUNDANT_COLUMNS),
    }


def shared_arguments(metadata):
    return {
        key: value
        for key, value in metadata.get("arguments", {}).items()
        if key not in SHARED_ARGUMENT_EXCLUSIONS
    }


def validate_shards(profile, shards):
    if not shards:
        raise ValueError("At least one shard is required")
    reference = shards[0]
    ref_meta = reference["metadata"]
    ref_args = shared_arguments(ref_meta)
    ref_grid_fields = reference["grid_fields"]
    ref_full_fields = reference["full_fields"]
    source_hash = ref_meta.get("source_sha256")
    all_grid_rows = []
    all_full_rows = []
    all_galaxies = set()
    shard_records = []

    for shard in shards:
        metadata = shard["metadata"]
        if metadata.get("profile") != profile:
            raise ValueError(f"Shard {shard['id']} profile does not match {profile!r}")
        if shared_arguments(metadata) != ref_args:
            raise ValueError(f"Shard {shard['id']} has different analysis settings")
        if metadata.get("source_sha256") != source_hash:
            raise ValueError(f"Shard {shard['id']} was produced by a different source version")
        if shard["grid_fields"] != ref_grid_fields:
            raise ValueError(f"Shard {shard['id']} profile-grid CSV schema differs")
        if shard["full_fields"] != ref_full_fields:
            raise ValueError(f"Shard {shard['id']} fully-marginalised CSV schema differs")

        shard_galaxies = {row["galaxy"] for row in shard["grid_rows"]}
        if not shard_galaxies:
            raise ValueError(f"Shard {shard['id']} contains no profile-grid galaxies")
        if len(shard_galaxies) != int(metadata.get("galaxy_count", -1)):
            raise ValueError(f"Shard {shard['id']} galaxy count disagrees with its metadata")
        overlap = all_galaxies.intersection(shard_galaxies)
        if overlap:
            raise ValueError(
                f"Galaxies occur in more than one shard: {', '.join(sorted(overlap))}"
            )
        all_galaxies.update(shard_galaxies)
        all_grid_rows.extend(shard["grid_rows"])
        all_full_rows.extend(shard["full_rows"])

        full_galaxies = {row["galaxy"] for row in shard["full_rows"]}
        if full_galaxies and full_galaxies != shard_galaxies:
            raise ValueError(f"Shard {shard['id']} has incomplete fully-marginalised rows")
        if not full_galaxies and metadata.get("disk_marginalisation", True):
            raise ValueError(f"Shard {shard['id']} is missing fully-marginalised results")

        shard_records.append({
            "shard_id": shard["id"],
            "galaxy_range": metadata.get("arguments", {}).get("galaxy_range"),
            "galaxy_count": len(shard_galaxies),
            "source_sha256": metadata.get("source_sha256"),
            "started_utc": metadata.get("started_utc"),
            "finished_utc": metadata.get("finished_utc"),
            "results_csv": str(shard["paths"]["profile_grid"]),
            "fully_marginalised_csv": str(shard["paths"]["fully_marginalised"]),
        })

    if all_full_rows and {row["galaxy"] for row in all_full_rows} != all_galaxies:
        raise ValueError("Combined fully-marginalised rows do not cover the profile-grid galaxies")
    return all_grid_rows, all_full_rows, ref_grid_fields, ref_full_fields, all_galaxies, shard_records


def write_csv_temp(directory, fieldnames, rows):
    with tempfile.NamedTemporaryFile(
        mode="w", newline="", encoding="utf-8", dir=directory,
        prefix=".merge-", suffix=".tmp", delete=False,
    ) as stream:
        temp_path = Path(stream.name)
        if fieldnames:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    return temp_path


def write_json_temp(directory, payload):
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=directory,
        prefix=".merge-", suffix=".tmp", delete=False,
    ) as stream:
        temp_path = Path(stream.name)
        json.dump(payload, stream, indent=2)
        stream.write("\n")
    return temp_path


def merge_shards(profile, output_root, shard_ids, replace_existing=False):
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    if len(shard_ids) != len(set(shard_ids)):
        raise ValueError("Shard IDs must be unique")
    if any(shard_id < 1 for shard_id in shard_ids):
        raise ValueError("Shard IDs must be positive integers")

    final_paths = {key: output_root / name for key, name in FINAL_NAMES.items()}
    existing = [path for path in final_paths.values() if path.exists()]
    if existing and not replace_existing:
        raise FileExistsError(
            "Merged outputs already exist; pass --replace-existing to replace: "
            + ", ".join(str(path) for path in existing)
        )

    shards = [read_shard(output_root, shard_id) for shard_id in sorted(shard_ids)]
    (grid_rows, full_rows, grid_fields, full_fields,
     galaxies, shard_records) = validate_shards(profile, shards)

    base = dict(shards[0]["metadata"])
    base["radius_coverage"] = {
        galaxy: coverage
        for shard in shards
        for galaxy, coverage in shard["metadata"].get("radius_coverage", {}).items()
    }

    base["arguments"] = shared_arguments(base)
    base["started_utc"] = min(
        shard["metadata"].get("started_utc", "") for shard in shards
    )
    base["finished_utc"] = max(
        shard["metadata"].get("finished_utc", "") for shard in shards
    )
    base["merged_utc"] = datetime.now(timezone.utc).isoformat()
    base["galaxy_count"] = len(galaxies)
    base["shard_count"] = len(shards)
    base["shard_ids"] = [shard["id"] for shard in shards]
    base["shards"] = shard_records
    base["results_csv"] = str(final_paths["profile_grid"])
    base["fully_marginalised_csv"] = str(final_paths["fully_marginalised"])

    temp_paths = []
    try:
        temp_paths.append(
            write_csv_temp(
                output_root,
                compact_fields(grid_fields, PROFILE_GRID_REDUNDANT_COLUMNS),
                grid_rows,
            )
        )
        temp_paths.append(
            write_csv_temp(
                output_root,
                compact_fields(full_fields, FULLY_MARGINALISED_REDUNDANT_COLUMNS),
                full_rows,
            )
        )
        temp_paths.append(write_json_temp(output_root, base))
        for temp_path, key in zip(temp_paths, FINAL_NAMES):
            os.replace(temp_path, final_paths[key])
    finally:
        for temp_path in temp_paths:
            temp_path.unlink(missing_ok=True)

    return len(galaxies), len(grid_rows), len(full_rows), final_paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("linear", "log_linear"), required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--shard-ids", type=int, nargs="+", required=True)
    parser.add_argument("--replace-existing", action="store_true")
    args = parser.parse_args()
    galaxies, grid_count, full_count, paths = merge_shards(
        args.profile, args.output_root, args.shard_ids, args.replace_existing)
    print(f"Merged {len(args.shard_ids)} shards covering {galaxies} galaxies")
    print(f"Profile-grid rows: {grid_count} -> {paths['profile_grid']}")
    print(f"Fully-marginalised rows: {full_count} -> {paths['fully_marginalised']}")
    print(f"Combined provenance: {paths['metadata']}")


if __name__ == "__main__":
    main()

"""Small standard-library tests for concurrent result-shard merging."""

import csv
import json
import tempfile
import unittest
from pathlib import Path

from merge_run_shards import merge_shards


class MergeRunShardsTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self):
        self.tempdir.cleanup()

    def write_shard(self, shard_id, galaxy, grid_size="17", profile="linear"):
        suffix = f".shard_{shard_id:02d}"
        grid_path = self.root / f"profile_grid{suffix}.csv"
        full_path = self.root / f"fully_marginalised{suffix}.csv"
        metadata_path = self.root / f"run_metadata{suffix}.json"
        rows = [{"galaxy": galaxy, "profile": profile, "log_probability": "-0.5"}]
        for path in (grid_path, full_path):
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=rows[0])
                writer.writeheader()
                writer.writerows(rows)
        metadata = {
            "profile": profile,
            "source_sha256": "same-source",
            "galaxy_count": 1,
            "disk_marginalisation": True,
            "started_utc": f"2026-09-23T10:0{shard_id}:00+00:00",
            "finished_utc": f"2026-09-23T10:1{shard_id}:00+00:00",
            "arguments": {
                "profile": profile,
                "grid_size": grid_size,
                "prob_grid_size": "128",
                "output_root": str(self.root),
                "galaxy_range": f"[{shard_id}, {shard_id}]",
                "shard_id": str(shard_id),
                "testing": "False",
            },
        }
        metadata_path.write_text(json.dumps(metadata))

    def test_merges_rows_and_records_each_shard(self):
        self.write_shard(1, "G1")
        self.write_shard(2, "G2")

        count, grid_count, full_count, paths = merge_shards(
            "linear", self.root, [1, 2])

        self.assertEqual((count, grid_count, full_count), (2, 2, 2))
        with paths["profile_grid"].open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.assertEqual({row["galaxy"] for row in rows}, {"G1", "G2"})
        metadata = json.loads(paths["metadata"].read_text())
        self.assertEqual(metadata["galaxy_count"], 2)
        self.assertEqual(metadata["shard_ids"], [1, 2])
        self.assertEqual(len(metadata["shards"]), 2)

    def test_missing_shard_is_rejected(self):
        self.write_shard(1, "G1")
        with self.assertRaises(FileNotFoundError):
            merge_shards("linear", self.root, [1, 2])

    def test_mismatched_settings_are_rejected(self):
        self.write_shard(1, "G1")
        self.write_shard(2, "G2", grid_size="21")
        with self.assertRaisesRegex(ValueError, "different analysis settings"):
            merge_shards("linear", self.root, [1, 2])

    def test_duplicate_galaxy_coverage_is_rejected(self):
        self.write_shard(1, "G1")
        self.write_shard(2, "G1")
        with self.assertRaisesRegex(ValueError, "more than one shard"):
            merge_shards("linear", self.root, [1, 2])

    def test_existing_aggregate_requires_explicit_replacement(self):
        self.write_shard(1, "G1")
        merge_shards("linear", self.root, [1])
        with self.assertRaises(FileExistsError):
            merge_shards("linear", self.root, [1])
        merge_shards("linear", self.root, [1], replace_existing=True)


if __name__ == "__main__":
    unittest.main()

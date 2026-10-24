"""Lightweight production GP interface tests; no Casertano integrals."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
import gp_rar_search as gp
import run_gp_rar_batch as batch


class BatchTests(unittest.TestCase):
    def test_exact_disjoint_seven_galaxy_allocation(self):
        items = [{"galaxy": f"G{i:03d}", "weight": (i+1)**2} for i in range(175)]
        shards = batch.allocate(items)
        self.assertEqual(len(shards), 25)
        self.assertTrue(all(len(s["galaxies"]) == 7 for s in shards))
        self.assertEqual(len({g for s in shards for g in s["galaxies"]}), 175)
        self.assertEqual(shards, batch.allocate(list(reversed(items))))

    def test_allocation_rejects_missing_and_duplicate_ids(self):
        with self.assertRaises(ValueError):
            batch.allocate([])
        with self.assertRaises(ValueError):
            batch.allocate([{"galaxy": "same", "weight": 1}]*175)

    def test_nested_direction_prefix_is_budget_independent(self):
        short = gp.nested_direction_rng(123, "NGC0100", .5, 0).standard_normal((5, 7))
        long = gp.nested_direction_rng(123, "NGC0100", .5, 0).standard_normal((10, 7))
        np.testing.assert_array_equal(short, long[:5])
        sequential = gp.nested_direction_rng(123, "NGC0100", .5, 0)
        np.testing.assert_array_equal(short, np.array([sequential.standard_normal(7) for _ in range(5)]))

    def test_streams_are_independent_of_galaxy_order_and_other_ell(self):
        a = gp.nested_direction_rng(123, "NGC0100", .5, 0).standard_normal(7)
        for galaxy, ell, replicate in [("NGC0100", 1., 0), ("D512-2", .5, 0), ("NGC0100", .5, 1)]:
            self.assertFalse(np.array_equal(a, gp.nested_direction_rng(123, galaxy, ell, replicate).standard_normal(7)))
        np.testing.assert_array_equal(a, gp.nested_direction_rng(123, "NGC0100", .5, 0).standard_normal(7))

    def test_compact_profile_roundtrip_and_collision_guard(self):
        profiles = [np.array([.5, .6]), np.array([.4, .7])]
        with tempfile.TemporaryDirectory() as folder:
            storage = gp.compact_profiles(Path(folder), profiles)
            array = np.load(storage["profiles_file"], mmap_mode="r", allow_pickle=False)
            np.testing.assert_array_equal(array, profiles)
            self.assertEqual(storage["sha256"], batch.sha(storage["profiles_file"]))
            with self.assertRaises(FileExistsError):
                gp.compact_profiles(Path(folder), profiles)

    def test_payload_is_explicit_and_queue_has_no_python_prefix(self):
        command = batch.payload("NGC0100", Path("/tmp/gp-test"))
        for flag, value in [("--directions", "5000"), ("--sampling-scheme", "nested-v1"),
                            ("--checkpoint-mode", "compact"), ("--ml-samples", "9")]:
            self.assertEqual(command[command.index(flag)+1], value)
        queue = batch.queue_commands(Path("/tmp/gp-test"))
        self.assertEqual(len(queue), 25)
        self.assertTrue(all('-q cmb -n 1 -m 8 -c "5 hrs"' in c for c in queue))
        self.assertTrue(all("python" not in c for c in queue))

    def test_shard_continues_after_process_failure(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            manifest = {"shards": [{"galaxies": ["bad", "good"]}], "references": {}}
            (output/"batch_manifest.json").write_text("{}")
            (output/"galaxies").mkdir()
            (output/"galaxies/good_gp_search.json").write_text(json.dumps({
                "status": "finite_search_completed", "records": [{"valid": True}]}))
            with patch.object(batch, "verify_manifest", return_value=manifest), patch.object(
                    batch.subprocess, "run", side_effect=[SimpleNamespace(returncode=1), SimpleNamespace(returncode=0)]) as run:
                with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
                    batch.run_shard(output, 1)
            self.assertEqual(run.call_count, 2)
            report = json.loads((output/"shard_01.json").read_text())
            self.assertEqual(report["failed_or_numerically_invalid"], 1)
            self.assertEqual(report["results"][1]["status"], "finite_search_completed")


    def test_cloud_reads_compact_nonfrontier_and_rejects_bad_hash(self):
        import plot_ml_profile_clouds as clouds
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            checkpoints = root/"G_gp_checkpoints"
            checkpoints.mkdir()
            storage = gp.compact_profiles(checkpoints, [np.array([.5, .5]), np.array([.4, .7])])
            data = {"status": "finite_search_completed", "galaxy": "G", "settings": {},
                    "provenance": {"source_sha256": {}, "script_sha256": "test", "input_sha256": {}},
                    "physical_domain": {}, "knots_kpc": [0., 1.], "frontier": [],
                    "checkpoint_storage": storage,
                    "records": [{"valid": True, "J_ref": 0., "marginal_logp": np.log(.01)},
                                {"valid": True, "J_ref": 3., "marginal_logp": np.log(.1)}]}
            path = root/"G_gp_search.json"
            path.write_text(json.dumps(data))
            panels, audit = clouds.load_gp([root], [2.], "one-sided")
            self.assertEqual(next(iter(panels.values()))[0]["ml"], [.4, .7])
            data["checkpoint_storage"]["sha256"] = "bad"
            path.write_text(json.dumps(data))
            with self.assertRaises(ValueError):
                clouds.load_gp([root], [2.], "one-sided")


if __name__ == "__main__":
    unittest.main()

"""Functional, synthetic-only tests for the versioned shard dispatcher."""
import hashlib, json, tempfile, unittest
from pathlib import Path
from unittest.mock import patch
import run_signed_rank_shard_20261019 as d

class DispatcherTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name)
        self.parent = self.base / "allowed"
        self.parent.mkdir()
        d.OUTPUT_PARENT = self.parent
        self.out = self.parent / "full_20261019_v1"
        self.runner = self.base / "runner_manifest.json"
        self.runner.write_text("{}")
        self.manifest_path = self.base / "dispatch.json"
        self.manifest = {"schema": d.SCHEMA, "signed_rank_policy": d.POLICY,
            "output_root": str(self.out), "runner_manifest": str(self.runner),
            "runner_manifest_sha256": hashlib.sha256(self.runner.read_bytes()).hexdigest(),
            "source_sha256": {str(Path(d.__file__).resolve()): d.sha256(d.__file__)},
            "shards": {"01": ["SYNTH-A"]}, "argv_template": ["fake", "--galaxy", "{galaxy}",
                "--output-root", "{output_root}", "--manifest", "{manifest}"],
            "max_runtime_seconds": 5, "minimum_free_bytes": 0}
        self.manifest_path.write_text(json.dumps(self.manifest))

    def tearDown(self):
        self.tmp.cleanup()

    def test_dry_run_does_not_create_output_or_run_science(self):
        commands = d.dispatch(self.manifest, self.manifest_path, "01", dry_run=True)
        self.assertEqual(commands[0][2], "SYNTH-A")
        self.assertFalse(self.out.exists())

    def test_duplicate_galaxy_allocation_rejected(self):
        self.manifest["shards"]["02"] = ["SYNTH-A"]
        with self.assertRaisesRegex(ValueError, "duplicates"):
            d.validate_manifest(self.manifest, self.manifest_path)

    def test_smoke_root_rejected(self):
        self.manifest["output_root"] = str(self.parent / "integration_smoke_v1_20261016")
        with self.assertRaisesRegex(ValueError, "pilot and production evidence roots are immutable"):
            d.validate_manifest(self.manifest, self.manifest_path)

    def test_success_claims_once_and_writes_status(self):
        class Result: returncode=0; stdout="synthetic"; stderr=""
        def fake_run(argv, **kwargs):
            output = Path(argv[argv.index("--output-root") + 1])
            output.mkdir(parents=True, exist_ok=True)
            (output / "SYNTH-A_gp_search.json").write_text("{}")
            return Result()
        with patch.object(d.subprocess, "run", side_effect=fake_run) as run:
            result = d.dispatch(self.manifest, self.manifest_path, "01", dry_run=False)
            self.assertEqual(result[0]["state"], "complete")
            run.assert_called_once()
        status = self.out / "SYNTH-A" / "dispatcher_status.json"
        self.assertEqual(json.loads(status.read_text())["state"], "complete")
        with self.assertRaises(FileExistsError):
            d.dispatch(self.manifest, self.manifest_path, "01", dry_run=False)

    def test_failed_runner_is_recorded_and_stops(self):
        class Result: returncode=3; stdout=""; stderr="synthetic failure"
        with patch.object(d.subprocess, "run", return_value=Result()):
            with self.assertRaisesRegex(RuntimeError, "without retry/refill"):
                d.dispatch(self.manifest, self.manifest_path, "01", dry_run=False)
        status = self.out / "SYNTH-A" / "dispatcher_status.json"
        self.assertEqual(json.loads(status.read_text())["state"], "failed")

    def reference_manifest(self):
        self.manifest.update(mode="reference", seeds={"SYNTH-A": 123},
            argv_template=["fake", "reference", "--galaxy", "{galaxy}",
                "--seed", "{seed}", "--output", "{output_file}"])

    def test_reference_seed_and_output_render_without_science(self):
        self.reference_manifest()
        command = d.dispatch(self.manifest, self.manifest_path, "01")[0]
        self.assertEqual(command[command.index("--seed")+1], "123")
        self.assertEqual(Path(command[-1]), self.out / "SYNTH-A" / "SYNTH-A_B1024.json")
        self.assertFalse(self.out.exists())

    def test_reference_requires_exact_seed_allocation(self):
        self.reference_manifest()
        self.manifest["seeds"] = {}
        with self.assertRaisesRegex(ValueError, "exact allocation"):
            d.dispatch(self.manifest, self.manifest_path, "01")

    def test_reference_success_uses_reference_result_filename(self):
        self.reference_manifest()
        class Result: returncode=0; stdout="synthetic reference"; stderr=""
        def fake_run(argv, **kwargs):
            Path(argv[-1]).write_text("{}")
            return Result()
        with patch.object(d.subprocess, "run", side_effect=fake_run) as run:
            result = d.dispatch(self.manifest, self.manifest_path, "01", dry_run=False)
        self.assertEqual(result[0]["state"], "complete")
        run.assert_called_once()

    def test_unknown_shard_does_not_create_output(self):
        with self.assertRaisesRegex(ValueError, "unknown shard"):
            d.dispatch(self.manifest, self.manifest_path, "99")
        self.assertFalse(self.out.exists())

if __name__ == "__main__": unittest.main()


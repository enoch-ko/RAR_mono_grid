"""Synthetic functional tests only: no Casertano/probability evaluations or queues."""
import contextlib, copy, importlib, io, json, os, sys, tempfile, types, unittest
from pathlib import Path
from unittest.mock import patch
import run_loglinear_bulgeless47_20261019 as r

class ReplacementTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.out = self.root/"science"
        self.audit = self.root/"audit.json"
        self.norm = self.root/"normalization.json"
        self.names = [f"G{i:03d}" for i in range(47)]
        counts = [10,10,13,14]
        start, jobs = 0, []
        for (job,sid),n in zip(r.JOBS.items(),counts):
            jobs.append(dict(job_id=job,shard=sid,unique_bulgeless=self.names[start:start+n]))
            start += n
        self.audit.write_text(json.dumps(dict(jobs=jobs)))
        self.norm.write_text(json.dumps(dict(galaxies={g:dict(Lbulge_1e9=0,kappa=1) for g in self.names})))
        self.m = dict(schema="loglinear-bulgeless47-v1",settings=dict(r.SETTINGS),output_root=str(self.out),
            normalization=str(self.norm),allocation_audit=str(self.audit),allocation_audit_sha256=r.sha(self.audit),
            shards={f"{i:02d}":[g] for i,g in enumerate(self.names,1)},verified_completed_galaxies=[],
            source_sha256={str(r.ENGINE):r.sha(r.ENGINE),str(Path(r.__file__).resolve()):r.sha(r.__file__)},
            input_sha256={str(self.norm):r.sha(self.norm)},data_root="/synthetic/data",dens_root="/synthetic/dens",table="/synthetic/table")
        self.output_patch = patch.object(r,"OUTPUT",self.out)
        self.output_patch.start()
    def tearDown(self):
        self.output_patch.stop()
        self.tmp.cleanup()
    def test_exact_allocation(self):
        self.assertEqual(r.validate(self.m),self.out)
        self.assertFalse(self.out.exists())
    def test_duplicate_or_missing_refused(self):
        self.m["shards"]["47"] = self.m["shards"]["46"]
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_multi_galaxy_shard_refused(self):
        self.m["shards"]["01"].append("G001")
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_other_namespace_refused(self):
        self.m["output_root"] = str(self.root/"other")
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_changed_settings_refused(self):
        self.m["settings"]["prob_grid_size"] = 800
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_input_pin_changed(self):
        self.norm.write_text("{}")
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_source_pin_changed(self):
        self.m["source_sha256"][str(r.ENGINE)] = "bad"
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_allocation_audit_pin_changed(self):
        self.audit.write_text("{}")
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_bulge_or_kappa_refused(self):
        for row in (dict(Lbulge_1e9=1,kappa=1),dict(Lbulge_1e9=0,kappa=.9)):
            d=json.loads(self.norm.read_text());d["galaxies"]["G000"]=row
            self.norm.write_text(json.dumps(d));self.m["input_sha256"][str(self.norm)]=r.sha(self.norm)
            with self.assertRaises(ValueError):r.validate(self.m)
    def test_unreviewed_exclusions_refused(self):
        self.m["verified_completed_galaxies"]=["G000"]
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_plan_flags_and_single_selection(self):
        with patch.object(r.engine,"MANIFEST",self.norm):
            p=r.build_plan(self.m,"01")
        args=p["argv"]
        self.assertEqual(args.count("--galaxy"),1)
        for key,value in r.SETTINGS.items():
            self.assertEqual(args[args.index("--"+key.replace("_","-"))+1],str(value))
        self.assertIn("--marginalise-gas-bulge",args)
        self.assertIn("--marginalise-disk",args)
        self.assertEqual(args[args.index("--shard-id")+1],"1")
        self.assertEqual(p["shard"],"01")
    def test_unknown_shard_refused(self):
        with self.assertRaises(ValueError):r.build_plan(self.m,"48")
    def test_claim_duplicate_and_durable_products(self):
        p=dict(output_root=str(self.out),shard="01",galaxies=["G000"])
        def mock_work(m,plan,event):
            for name,suffix in (("profile_grid","csv"),("fully_marginalised","csv"),("run_metadata","json")):
                (self.out/(name+".shard_01."+suffix)).write_text("synthetic fixture")
        r.run(self.m,p,"synthetic-job",work=mock_work)
        self.assertTrue((self.out/"claims/shard_01.json").exists())
        self.assertIn("shard_complete",(self.out/"shard_01.events.jsonl").read_text())
        with self.assertRaises(FileExistsError):r.run(self.m,p,"second",work=mock_work)
    def test_existing_claim_blocks_no_products_retry(self):
        p=dict(output_root=str(self.out),shard="01",galaxies=["G000"])
        def failed(*args):raise RuntimeError("fixture")
        with self.assertRaises(RuntimeError):r.run(self.m,p,"first",work=failed)
        with self.assertRaises(FileExistsError):r.run(self.m,p,"second",work=failed)
        self.assertIn("shard_failed",(self.out/"shard_01.events.jsonl").read_text())
    def test_missing_products_not_complete(self):
        with self.assertRaises(RuntimeError):
            r.run(self.m,dict(output_root=str(self.out),shard="01"),"fixture",work=lambda *a:None)
    def test_scheduler_guard_no_output(self):
        with self.assertRaises(ValueError):r.run(self.m,dict(output_root=str(self.out)),None)
        self.assertFalse(self.out.exists())
    def test_private_engine_does_not_mutate_public_module(self):
        public=importlib.import_module("run_loglinear_recovery_20261019")
        self.assertNotEqual(public.MANIFEST,r.engine.MANIFEST)
        self.assertIsNot(public,r.engine)
        self.assertEqual(public.PREP.name,"preparation_v1")
    def test_worker_mock_install_once_restores_hooks(self):
        self.out.mkdir()
        meta=self.out/"run_metadata.shard_01.json"
        mono=types.ModuleType("mono_SPARC_RAR")
        original=lambda galaxy,*a,**kw: "synthetic-result"
        mono.process_galaxy=original
        boundary=types.ModuleType("mono_SPARC_RAR_boundary")
        original_hashes=lambda:{"source":"fixture"}
        boundary.hashes=original_hashes
        events=[]
        def fake_main():
            self.assertEqual(mono.process_galaxy("G000"),"synthetic-result")
            meta.write_text("{}")
        boundary.main=fake_main
        recovery=types.ModuleType("run_corrected_recovery_pilot")
        calls=[]
        recovery.normalization_proof=lambda m,p:{"fixture":True}
        recovery.install=lambda m,p:calls.append(p)
        with patch.dict(sys.modules,{"mono_SPARC_RAR":mono,"mono_SPARC_RAR_boundary":boundary,
                "run_corrected_recovery_pilot":recovery}),patch.object(r.engine,"load_context",lambda:self.m),patch.object(r.engine,"MANIFEST",self.norm):
            r.engine.worker(self.m,dict(argv=["fixture"],output_root=str(self.out),shard="01",galaxies=["G000"]),events.append)
        self.assertEqual(len(calls),1)
        self.assertIs(mono.process_galaxy,original)
        self.assertIs(boundary.hashes,original_hashes)
        self.assertEqual([e["state"] for e in events],["galaxy_start","galaxy_finished"])
        self.assertTrue(json.loads(meta.read_text())["recovery_scope"]["normalization_applied_once"])

if __name__=="__main__":unittest.main()

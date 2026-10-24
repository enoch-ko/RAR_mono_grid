"""Synthetic/mocked checks only; no real-galaxy force, probability or grid runs."""
import contextlib, copy, io, json, os, sys, tempfile, types, unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
import run_loglinear_recovery_20261019 as d

class RecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.base=Path(self.temp.name)
        self.parent=self.base/"allowed";self.parent.mkdir()
        self.manifest=self.base/"manifest.json";self.manifest.write_text("{}")
        records={"B"+f"{i:02d}":dict(Ltotal_1e9=1.,Lbulge_1e9=.2,Ldisk_1e9=.8,kappa=.8) for i in range(32)}
        records.update({"C"+f"{i:02d}":dict(Ltotal_1e9=1.,Lbulge_1e9=0.,Ldisk_1e9=1.,kappa=1.) for i in range(27)})
        self.norm=self.base/"normalization.json"
        self.norm.write_text(json.dumps(dict(policy="disc_to_catalogue_total_minus_published_bulge_v1",galaxies=records)))
        names=sorted(records);groups={};offset=0
        for i in range(22):
            size=3 if i<15 else 2;groups[f"{i+1:02d}"]=names[offset:offset+size];offset+=size
        self.m=dict(schema="loglinear-recovery-v1",settings=dict(d.SETTINGS),normalization=str(self.norm),
            output_root=str(self.parent/"log_linear_v1"),failed_bulgeless_galaxies=[g for g in names if g.startswith("C")],
            shards=groups,source_sha256={str(Path(d.__file__).resolve()):d.sha(d.__file__)},
            input_sha256={str(self.norm):d.sha(self.norm)},data_root="/synthetic/data",dens_root="/synthetic/dens",table="/synthetic/table")
        self.patches=[patch.object(d,"OUTPUT_PARENT",self.parent),patch.object(d,"MANIFEST",self.manifest)]
        for p in self.patches:p.start()
    def tearDown(self):
        for p in reversed(self.patches):p.stop()
        self.temp.cleanup()
    def test_exact_allocation_and_no_output_dry_plan(self):
        d.validate(self.m);plan=d.build_plan(self.m,"01")
        self.assertEqual(plan["galaxies"],["B00","B01","B02"])
        self.assertFalse(Path(self.m["output_root"]).exists())
        self.assertIn("--marginalise-disk",plan["argv"])
        self.assertIn("--marginalise-gas-bulge",plan["argv"])
    def test_real_parser_accepts_all_planned_flags_without_calculation(self):
        import mono_SPARC_RAR as mono
        plan=d.build_plan(self.m,"01")
        # Boundary wrapper consumes its two refinement options before core parsing.
        argv=list(plan["argv"])
        for flag in ("--boundary-log-tolerance","--boundary-max-refinements"):
            j=argv.index(flag);del argv[j:j+2]
        with patch.object(sys,"argv",argv):args=mono.parse_args()
        self.assertEqual((args.profile,args.grid_size,args.prob_grid_size,args.ml_samples),("log_linear",51,400,9))
        self.assertEqual(args.galaxy,plan["galaxies"])
        self.assertEqual(args.gas_ml,1.);self.assertEqual(args.bulge_ml,.7)
    def test_duplicate_allocation_rejected(self):
        self.m["shards"]["02"][0]="B00"
        with self.assertRaisesRegex(ValueError,"unique galaxies"):d.validate(self.m)
    def test_changed_source_pin_rejected(self):
        self.m["source_sha256"][str(Path(d.__file__).resolve())]="0"*64
        with self.assertRaisesRegex(ValueError,"dependency changed"):d.validate(self.m)
    def test_login_execution_rejected_before_context_or_outputs(self):
        with patch.dict(os.environ,{},clear=True),patch.object(sys,"argv",["worker","01"]),patch.object(d,"load_context") as load:
            with contextlib.redirect_stderr(io.StringIO()),self.assertRaises(SystemExit):d.main()
            load.assert_not_called()
        self.assertFalse(Path(self.m["output_root"]).exists())
    def test_mock_success_is_exclusive_and_records_failure(self):
        plan=d.build_plan(self.m,"01")
        def fake(m,plan,event):
            for prefix,suffix in (("profile_grid","csv"),("fully_marginalised","csv"),("run_metadata","json")):
                (Path(plan["output_root"])/(prefix+".shard_01."+suffix)).write_text("{}")
        with contextlib.redirect_stdout(io.StringIO()):d.run(self.m,plan,"SYNTHETIC",work=fake)
        with self.assertRaises(FileExistsError):d.run(self.m,plan,"SYNTHETIC",work=fake)
        plan=d.build_plan(self.m,"02")
        def failing(*unused):raise RuntimeError("synthetic worker failure")
        with contextlib.redirect_stdout(io.StringIO()),self.assertRaisesRegex(RuntimeError,"synthetic"):
            d.run(self.m,plan,"SYNTHETIC",work=failing)
        events=[json.loads(s) for s in (Path(plan["output_root"])/"shard_02.events.jsonl").read_text().splitlines()]
        self.assertEqual(events[-1]["state"],"shard_failed")
    def test_actual_normalization_adapter_called_once_with_mocked_boundary(self):
        rad=np.array([1.,2.]);raw=np.ones(2)*(1e9/(3*np.pi));seen=[]
        utils=types.ModuleType("utils");utils.DEFAULT_DENS_ROOT=Path("/synthetic/dens")
        utils.get_SBdisk=lambda *a,**kw:(rad.copy(),raw.copy())
        utils.load_sparc_table=lambda:pd.DataFrame({"Galaxy":["B00","C00"],"L":[1.,1.]})
        utils.compute_A_B=lambda r,s:(2*np.pi*np.trapz(s*r,r),0.)
        mono=types.ModuleType("mono_SPARC_RAR");mono.get_SBdisk=utils.get_SBdisk
        inverse=types.ModuleType("inverse_rar_diagnostic_corrected");inverse.get_SBdisk=utils.get_SBdisk;inverse.cache_gate=lambda *a:{}
        def process(g,*unused,**kw):
            seen.append(mono.get_SBdisk(g)[1]);seen.append(mono.get_SBdisk(g)[1])
        mono.process_galaxy=process
        boundary=types.ModuleType("mono_SPARC_RAR_boundary");boundary.hashes=lambda:{}
        plan=d.build_plan(self.m,"01");root=Path(plan["output_root"]);root.mkdir()
        def main():
            mono.process_galaxy("B00");mono.process_galaxy("C00")
            (root/"run_metadata.shard_01.json").write_text("{}")
        boundary.main=main
        modules={"utils":utils,"mono_SPARC_RAR":mono,"inverse_rar_diagnostic_corrected":inverse,"mono_SPARC_RAR_boundary":boundary}
        with patch.dict(sys.modules,modules),patch.object(d,"load_context",return_value=self.m):
            d.worker(self.m,plan,lambda event:None)
        np.testing.assert_allclose(seen[0],raw*.8);np.testing.assert_array_equal(seen[0],seen[1])
        np.testing.assert_allclose(seen[2],raw);np.testing.assert_array_equal(seen[2],seen[3])
        metadata=json.loads((root/"run_metadata.shard_01.json").read_text())
        self.assertTrue(metadata["recovery_scope"]["normalization_applied_once"])
if __name__=="__main__":unittest.main()

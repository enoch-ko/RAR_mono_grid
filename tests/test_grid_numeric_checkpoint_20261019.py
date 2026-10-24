"""Synthetic checkpoint/launcher tests; no real force, probability, or queue calls."""
import copy,csv,importlib,json,math,sys,tempfile,types,unittest
from pathlib import Path
from unittest.mock import patch
import grid_numeric_checkpoint_20261019 as c
import run_loglinear_checkpoint37_20261019 as r

def axes():return dict(log10_mean_ml=[-1.,0.],mean_ml=[.1,1.],slopes=[-.2,.2])
def pair(i,j,galaxy="G000",value=-12000.25,valid=True):
    a=axes()
    return [dict(galaxy=galaxy,profile="log_linear",probability_kind=kind,grid_i=i,grid_j=j,
        log10_upsilon_disk=a["log10_mean_ml"][j],upsilon_disk_mean=a["mean_ml"][j],
        disk_slope=a["slopes"][i],log_probability=value,valid=valid,
        invalid_reason="" if valid else "fixture_invalid",warnings="exact warning | fixture")
        for kind in ("conditional","marginalised_gas_bulge")]
def read(path):
    with Path(path).open() as stream:return list(csv.DictReader(stream))

class CheckpointTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.identity=dict(galaxy="G000",manifest_sha256="fixture",source_sha256={"x":"y"})
        self.cp=c.GridCheckpoint(self.root/"checkpoint",self.identity,axes())
    def tearDown(self):
        for stream in self.cp.streams.values():
            if not stream.closed:stream.close()
        self.tmp.cleanup()
    def fill(self):
        for j in range(2):
            for i in range(2):self.cp.append_cell(pair(i,j))
    def finish(self):
        self.fill();self.cp.finish_grids([[False,False],[False,False]])
        self.cp.save_fully(dict(galaxy="G000",log_probability=-12000.25,valid=True))
        self.cp.finish_galaxy(dict(retained=4))
    def test_stream_flushed_before_grid_end(self):
        self.cp.append_cell(pair(0,0))
        rows=read(self.cp.path/"conditional.partial.csv")
        self.assertEqual(rows[0]["log_probability"],"-12000.25")
        self.assertFalse((self.cp.path/"conditional.csv").exists())
        self.assertFalse((self.cp.path/"grids_complete.json").exists())
    def test_complete_exact_values_hashes(self):
        self.finish();report=c.verify_complete(self.cp.path,self.identity)
        self.assertEqual(report["cells"],4)
        self.assertFalse(report["resume_implemented"])
        self.assertFalse(report["science_validated"])
        self.assertEqual(len(read(self.cp.path/"profile_grid.csv")),8)
        self.assertEqual([row["probability_kind"] for row in read(self.cp.path/"profile_grid.csv")],
            ["conditional","marginalised_gas_bulge"]*4)
    def test_no_clipping_invalid_reason_and_nan_retained(self):
        self.cp.append_cell(pair(0,0,value=float("nan"),valid=False))
        row=read(self.cp.path/"conditional.partial.csv")[0]
        self.assertEqual(row["log_probability"],"nan")
        self.assertEqual(row["valid"],"False")
        self.assertEqual(row["invalid_reason"],"fixture_invalid")
        self.assertEqual(row["warnings"],"exact warning | fixture")
    def test_duplicate_refused(self):
        self.cp.append_cell(pair(0,0))
        with self.assertRaises(ValueError):self.cp.append_cell(pair(0,0))
    def test_incomplete_grid_not_published(self):
        self.cp.append_cell(pair(0,0))
        with self.assertRaises(ValueError):self.cp.finish_grids([[False]*2]*2)
        self.assertFalse((self.cp.path/"grids_complete.json").exists())
    def test_mismatched_pair_refused(self):
        rows=pair(0,0);rows[1]["grid_j"]=1
        with self.assertRaises(ValueError):self.cp.append_cell(rows)
    def test_wrong_axes_refused(self):
        rows=pair(0,0);rows[0]["disk_slope"]=99
        with self.assertRaises(ValueError):self.cp.append_cell(rows)
    def test_wrong_galaxy_refused(self):
        with self.assertRaises(ValueError):self.cp.append_cell(pair(0,0,galaxy="OTHER"))
    def test_failed_partial_preserved(self):
        self.cp.append_cell(pair(0,0));self.cp.fail(RuntimeError("fixture interruption"))
        self.assertEqual(len(read(self.cp.path/"conditional.partial.csv")),1)
        failure=json.loads((self.cp.path/"failure.json").read_text())
        self.assertEqual(failure["cells_flushed"],1)
        self.assertFalse(failure["fully_marginalised_saved"])
        with self.assertRaises(FileNotFoundError):c.verify_complete(self.cp.path,self.identity)
    def test_grids_saved_before_full(self):
        self.fill();self.cp.finish_grids([[False]*2]*2)
        self.assertEqual(len(read(self.cp.path/"conditional.csv")),4)
        self.assertFalse((self.cp.path/"fully_marginalised.csv").exists())
        self.assertFalse((self.cp.path/"galaxy_complete.json").exists())
        self.cp.fail(RuntimeError("full integral fixture failure"))
        self.assertTrue(json.loads((self.cp.path/"failure.json").read_text())["grids_ready"])
    def test_partial_atomic_publication_failure_not_complete(self):
        self.fill()
        with patch.object(c,"publish_json",side_effect=OSError("fixture quota")):
            with self.assertRaises(OSError):self.cp.finish_grids([[False]*2]*2)
        self.assertFalse((self.cp.path/"grids_complete.json").exists())
        self.assertFalse((self.cp.path/"galaxy_complete.json").exists())
    def test_no_overwrite_or_existing_directory(self):
        with self.assertRaises(FileExistsError):c.GridCheckpoint(self.cp.path,self.identity,axes())
        with self.assertRaises(FileExistsError):c.publish_json(self.cp.path/"context.json",{})
    def test_corruption_rejected(self):
        self.finish();(self.cp.path/"conditional.csv").write_text("corrupted fixture")
        with self.assertRaises(ValueError):c.verify_complete(self.cp.path,self.identity)
    def test_identity_change_rejected(self):
        self.finish()
        with self.assertRaises(ValueError):c.verify_complete(self.cp.path,dict(self.identity,manifest_sha256="other"))
    def test_full_cannot_precede_grid(self):
        with self.assertRaises(ValueError):self.cp.save_fully({"valid":True})
        with self.assertRaises(ValueError):self.cp.finish_galaxy({})
    def test_json_nonfinite_policy(self):
        self.assertIsNone(c.json_value({"p":float("nan")})["p"])
    def test_published_alias_same_inode(self):
        self.finish()
        self.assertEqual((self.cp.path/"conditional.csv").stat().st_ino,
            (self.cp.path/"conditional.partial.csv").stat().st_ino)

class ScopeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name);self.out=self.root/"output"
        self.audit=self.root/"audit.json";self.norm=self.root/"norm.json"
        names=[f"G{i:03d}" for i in range(47)]
        self.audit.write_text(json.dumps(dict(jobs=[dict(job_id=j,shard=s,unique_bulgeless=group)
            for j,s,group in ((935526,5,names[:10]),(935528,7,names[10:20]),
                (935529,8,names[20:33]),(935531,10,names[33:]))])))
        self.norm.write_text(json.dumps(dict(galaxies={g:dict(Lbulge_1e9=0,kappa=1) for g in names})))
        self.m=dict(schema="loglinear-checkpoint37-v1",settings=dict(r.SETTINGS),output_root=str(self.out),
            allocation_audit=str(self.audit),allocation_audit_sha256=r.sha(self.audit),normalization=str(self.norm),
            shards={f"{i:02d}":[g] for i,g in enumerate(names[10:],1)},keep_job=935526,
            replacing_jobs=list(r.JOBS),verified_completed_galaxies=[],
            source_sha256={str(p.resolve()):r.sha(p) for p in (r.ENGINE,Path(r.__file__),Path(c.__file__))},
            input_sha256={str(self.norm):r.sha(self.norm)},data_root="/fixture/data",dens_root="/fixture/dens",table="/fixture/table")
        self.outpatch=patch.object(r,"OUTPUT",self.out);self.outpatch.start()
    def tearDown(self):self.outpatch.stop();self.tmp.cleanup()
    def test37_excludes_all10_kept(self):self.assertEqual(r.validate(self.m),self.out)
    def test_kept_galaxy_forbidden(self):
        self.m["shards"]["01"]=["G000"]
        with self.assertRaises(ValueError):r.validate(self.m)
    def test47_old_scope_refused(self):
        self.m["schema"]="loglinear-bulgeless47-v1"
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_duplicate_target_refused(self):
        self.m["shards"]["01"]=self.m["shards"]["02"]
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_changed_pin_or_settings(self):
        self.m["settings"]["grid_size"]=31
        with self.assertRaises(ValueError):r.validate(self.m)
        self.m["settings"]=dict(r.SETTINGS);self.m["source_sha256"][str(r.ENGINE)]="bad"
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_changed_keep_or_cancel(self):
        self.m["replacing_jobs"].append(935526)
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_changed_checkpoint_source(self):
        self.m["source_sha256"][str(Path(c.__file__).resolve())]="bad"
        with self.assertRaises(ValueError):r.validate(self.m)
    def test_plan_single_galaxy_and_flags(self):
        with patch.object(r.engine,"MANIFEST",self.norm):plan=r.build_plan(self.m,"01")
        self.assertEqual(plan["galaxies"],["G010"])
        self.assertEqual(plan["argv"].count("--galaxy"),1)
        for key,value in r.SETTINGS.items():
            self.assertEqual(plan["argv"][plan["argv"].index("--"+key.replace("_","-"))+1],str(value))
        self.assertIn("--marginalise-disk",plan["argv"])
    def test_unknown_shard_and_login_refused(self):
        with self.assertRaises(ValueError):r.build_plan(self.m,"38")
        with self.assertRaises(ValueError):r.run(self.m,dict(output_root=str(self.out)),None)
        self.assertFalse(self.out.exists())
    def test_exclusive_claim_and_failure_persists(self):
        plan=dict(output_root=str(self.out),shard="01",galaxies=["G010"])
        def failed(*args):raise RuntimeError("fixture")
        with self.assertRaises(RuntimeError):r.run(self.m,plan,"fixture-job",work=failed)
        with self.assertRaises(FileExistsError):r.run(self.m,plan,"second",work=failed)
        self.assertIn("shard_failed",(self.out/"shard_01.events.jsonl").read_text())
    def test_frozen_public_worker_unchanged(self):
        public=importlib.import_module("run_loglinear_recovery_20261019")
        self.assertIsNot(public,r.engine);self.assertEqual(public.PREP.name,"preparation_v1")
        self.assertEqual(public.validate.__module__,public.__name__)

class HookTests(unittest.TestCase):
    def exercise(self,fail_after=None):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/"out";mf=Path(tmp)/"manifest.json";mf.write_text("{}")
            m=dict(source_sha256={"fixture":"source"},input_sha256={},settings=dict(r.SETTINGS),normalization=str(mf))
            plan=dict(output_root=str(out),shard="01",galaxies=["G000"],argv=["fixture"])
            mono=types.ModuleType("mono_SPARC_RAR");boundary=types.ModuleType("mono_SPARC_RAR_boundary")
            recovery=types.ModuleType("run_corrected_recovery_pilot")
            calls=dict(install=0,cell=0,full=0);events=[];saved_rows=[]
            def install(*args):calls["install"]+=1
            recovery.install=install;recovery.normalization_proof=lambda *args:{"fixture":True}
            def append(rows,galaxy,profile,i,j,*args):rows.extend(pair(i,j))
            def full(*args):
                calls["full"]+=1
                self.assertTrue((out/"checkpoints/G000/grids_complete.json").exists())
                return dict(galaxy="G000",log_probability=-12000.25,valid=True)
            def process(galaxy,table,args,rows,x,means,slopes,*nuisance):
                for j in range(2):
                    for i in range(2):
                        mono.append_probability_rows(rows,galaxy,"log_linear",i,j)
                        calls["cell"]+=1
                        if fail_after==calls["cell"]:raise RuntimeError("fixture interruption")
                saved_rows.extend(rows)
                f=mono.fully_marginalise_grid(galaxy,args,x,slopes,[],[],[[False]*2]*2)
                return [],[],[],f,dict(retained=4)
            mono.append_probability_rows=append;mono.fully_marginalise_grid=full;mono.process_galaxy=process
            original_hashes=lambda:{"fixture":"hash"};boundary.hashes=original_hashes
            def main():
                a=axes()
                mono.process_galaxy("G000",None,types.SimpleNamespace(),[],a["log10_mean_ml"],a["mean_ml"],a["slopes"],[1],[1],[.7],[1])
                (out/"run_metadata.shard_01.json").write_text('{"preserved":true}')
            boundary.main=main
            with patch.dict(sys.modules,{"mono_SPARC_RAR":mono,"mono_SPARC_RAR_boundary":boundary,
                    "run_corrected_recovery_pilot":recovery}),patch.object(r,"load_context",lambda:m),patch.object(r,"MANIFEST",mf):
                if fail_after is None:r.worker(m,plan,events.append)
                else:
                    with self.assertRaisesRegex(RuntimeError,"fixture interruption"):r.worker(m,plan,events.append)
            self.assertIs(mono.process_galaxy,process);self.assertIs(mono.append_probability_rows,append)
            self.assertIs(mono.fully_marginalise_grid,full);self.assertIs(boundary.hashes,original_hashes)
            self.assertEqual(calls["install"],1)
            cp=out/"checkpoints/G000"
            if fail_after is None:
                context=json.loads((cp/"context.json").read_text());c.verify_complete(cp,context["identity"])
                self.assertEqual(calls["cell"],4);self.assertEqual(calls["full"],1)
                actual=read(cp/"profile_grid.csv")
                self.assertEqual(actual,[{key:str(value) for key,value in row.items()} for row in saved_rows])
                self.assertEqual(json.loads((out/"run_metadata.shard_01.json").read_text()),{"preserved":True})
                self.assertTrue((out/"checkpoint_provenance.shard_01.json").exists())
                self.assertEqual([e["state"] for e in events],["galaxy_start","numeric_grids_checkpointed","galaxy_finished"])
            else:
                self.assertEqual(len(read(cp/"conditional.partial.csv")),fail_after)
                self.assertFalse((cp/"galaxy_complete.json").exists())
                self.assertEqual(calls["full"],0)
    def test_mock_original_cell_order_and_values_unchanged(self):self.exercise()
    def test_mock_failure_saves_completed_cells_restores_hooks(self):self.exercise(fail_after=2)

if __name__=="__main__":unittest.main()

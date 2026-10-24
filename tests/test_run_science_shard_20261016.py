"""Packaging-only tests. Scientific child execution is MOCKED; addqueue never runs."""
import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock,patch
import run_science_shard_20261016 as worker


class PackagingTests(unittest.TestCase):
    def temporary_plan(self,root):
        plan=worker.build_plan('gp',1,['G'+str(i) for i in range(7)])
        plan['receipt_root']=str(root/'receipts');plan['log_root']=str(root/'logs')
        for item in plan['entries']:
            item['output']=str(root/(item['galaxy']+'.json'))
            item['collision_guards']=[item['output']]
        return plan

    def run_mock(self,plan,runner,clock=worker.time.monotonic):
        with contextlib.redirect_stdout(io.StringIO()):
            return worker.run_plan(plan,'MOCK_JOB',runner=runner,clock=clock)

    def test_python_syntax(self):
        compile(Path(worker.__file__).read_text(),str(worker.__file__),'exec')

    def test_config_and_manifest(self):
        config,manifest,allocation=worker.load_context()
        self.assertEqual(config['worker_sha256'],worker.sha(worker.__file__))
        self.assertEqual(worker.sha(worker.MANIFEST),worker.MANIFEST_SHA)
        self.assertEqual(len(allocation),25)

    def test_exact_25_scope(self):
        self.assertEqual(sum(map(len,worker.ALLOWED.values())),25)
        self.assertEqual(worker.ALLOWED['gp'],list(range(1,21)))
        self.assertEqual(worker.ALLOWED['reference'],[1,6,8,9,15])

    def test_out_of_scope_and_injection(self):
        for mode,shard in [('gp','21'),('reference','2'),('qp','1'),('gp','1;touch'),('gp','-1')]:
            with self.assertRaises(ValueError):worker.mode_shard(mode,shard)

    def test_no_duplicate_targets(self):
        _,_,allocation=worker.load_context()
        by_mode={}
        for mode,ids in worker.ALLOWED.items():
            entries=[e for shard in ids for e in worker.build_plan(mode,shard,allocation[shard])['entries']]
            targets=[e['output'] for e in entries]
            self.assertEqual(len(targets),len(set(targets)))
            by_mode[mode]=set(targets)
        self.assertTrue(by_mode['gp'].isdisjoint(by_mode['reference']))
        self.assertEqual(len(by_mode['gp']),140)
        self.assertEqual(len(by_mode['reference']),35)

    def test_seeds(self):
        _,_,allocation=worker.load_context()
        seeds=[e['seed'] for s in worker.ALLOWED['reference'] for e in worker.build_plan('reference',s,allocation[s])['entries']]
        expected=[202610160000+7*(s-1)+i for s in [1,6,8,9,15] for i in range(1,8)]
        self.assertEqual(seeds,expected)
        self.assertEqual(len(set(seeds)),35)

    def test_simple_queue_serialization(self):
        for mode,ids in worker.ALLOWED.items():
            for shard in ids:
                tokens=[str(Path(worker.__file__).resolve()),mode,format(shard,'02d')]
                self.assertEqual(shlex.split(' '.join(tokens)),tokens)
                check=subprocess.run(['/bin/bash','-n','-c',' '.join(tokens)],capture_output=True)
                self.assertEqual(check.returncode,0)

    def test_old_failure_syntax_only(self):
        old=worker.SCIENCE/'logs/timeout.sh'
        result=subprocess.run(['/bin/bash','-n',str(old)],capture_output=True,text=True)
        self.assertEqual(result.returncode,2)
        self.assertIn('unexpected token',result.stderr)

    def test_gp_exact_arguments(self):
        args=worker.build_argv('gp','UGC07399')
        expected=[str(worker.ROOT/'run_corrected_recovery_pilot.py'),'gp','--manifest',str(worker.MANIFEST),
            '--seconds','172740','--','--galaxy','UGC07399','--output-root',str(worker.SCIENCE/'gp_d5000_a16'),
            '--data-root','/mnt/users/koe/data','--dens-root','/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG',
            '--validation-root',str(worker.ROOT/'outputs'),'--directions','5000','--amplitude-steps','16',
            '--ell-over-rd','0.5','1','2','4','--length-fraction','0.3','--min-ml','0.01','--max-ml','10',
            '--gbar-floor-relative','1e-8','--seed','20261016','--sampling-scheme','nested-v1','--replicate','0',
            '--checkpoint-mode','retained-only','--max-retained-bytes','8388608',
            '--retain-probability','0.0013498980316300957','--gas-ml','1','--gas-ml-dex','0.04',
            '--bulge-ml','0.7','--bulge-ml-dex','0.1','--ml-samples','9','--prob-grid-size','400',
            '--zpadding','8','--boundary-log-tolerance','2e-5','--boundary-max-refinements','8']
        self.assertEqual(args,expected)
        self.assertNotIn('--reference-json',args)

    def test_reference_exact_arguments_and_filename(self):
        args=worker.build_argv('reference','NGC0891',202610160038)
        expected=[str(worker.ROOT/'run_corrected_recovery_pilot.py'),'reference','--manifest',str(worker.MANIFEST),
            '--seconds','21540','--','build','--galaxy','NGC0891','--B','1024','--seed','202610160038',
            '--output',str(worker.SCIENCE/'reference_B1024/NGC0891_B1024.json'),
            '--normalization',str(worker.ROOT/'outputs/bulge_normalized_20261015/inputs/normalization.json'),
            '--data-root','/mnt/users/koe/data','--dens-root','/mnt/users/koe/SPARC_RAR/BulgeDiskDec_LTG',
            '--table',str(worker.ROOT/'SPARC_table.txt'),'--prob-grid-size','400','--zpadding','8',
            '--log-tolerance','2e-5','--max-refinements','8','--seconds','21540','--save-vectors']
        self.assertEqual(args,expected)
        self.assertNotIn('$',' '.join(args))

    def test_caps_and_resources(self):
        self.assertEqual(worker.LIMITS,{'gp':{'per_galaxy':172800,'wrapper':172740,'shard':1209600,'memory_gib':8},
            'reference':{'per_galaxy':21600,'wrapper':21540,'shard':172800,'memory_gib':4}})

    def test_success_receipts(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);plan=self.temporary_plan(root);calls=[]
            def fake(argv,seconds,log):
                calls.append((argv,seconds,log));output=next(x['output'] for x in plan['entries'] if x['argv']==argv)
                Path(output).write_text('{"MOCK":true}')
                return {'returncode':0,'timed_out':False,'elapsed_seconds':0.}
            self.assertEqual(self.run_mock(plan,fake),0)
            self.assertEqual(len(calls),7)
            summary=json.loads((root/'receipts/gp_01.summary.json').read_text())
            self.assertEqual(len(summary['results']),7)
            self.assertTrue(all(x['status']=='completed' for x in summary['results']))

    def test_failure_continues_other_galaxies(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t);plan=self.temporary_plan(root);calls=[]
            def fake(argv,seconds,log):
                calls.append(argv)
                if len(calls)>1:
                    Path(next(x['output'] for x in plan['entries'] if x['argv']==argv)).write_text('{"MOCK":true}')
                return {'returncode':1 if len(calls)==1 else 0,'timed_out':False,'elapsed_seconds':0.}
            self.assertEqual(self.run_mock(plan,fake),1)
            self.assertEqual(len(calls),7)

    def test_missing_output_not_success(self):
        with tempfile.TemporaryDirectory() as t:
            plan=self.temporary_plan(Path(t))
            fake=Mock(return_value={'returncode':0,'timed_out':False,'elapsed_seconds':0.})
            self.assertEqual(self.run_mock(plan,fake),1)
            self.assertEqual(fake.call_count,7)

    def test_science_collision_stops(self):
        with tempfile.TemporaryDirectory() as t:
            plan=self.temporary_plan(Path(t));Path(plan['entries'][0]['output']).write_text('MOCK')
            fake=Mock()
            with self.assertRaises(FileExistsError):self.run_mock(plan,fake)
            fake.assert_not_called()

    def test_claim_stops_duplicate_attempt(self):
        with tempfile.TemporaryDirectory() as t:
            plan=self.temporary_plan(Path(t));fake=Mock(return_value={'returncode':1,'timed_out':False})
            self.run_mock(plan,fake);fake.reset_mock()
            with self.assertRaises(FileExistsError):self.run_mock(plan,fake)
            fake.assert_not_called()

    def test_io_error_stops_instead_of_continuing(self):
        with tempfile.TemporaryDirectory() as t:
            plan=self.temporary_plan(Path(t));fake=Mock(side_effect=OSError('MOCK quota error'))
            with self.assertRaises(OSError):self.run_mock(plan,fake)
            self.assertEqual(fake.call_count,1)

    def test_shard_cap_no_child(self):
        with tempfile.TemporaryDirectory() as t:
            plan=self.temporary_plan(Path(t));fake=Mock();values=iter([0]+[1209601]*8)
            self.assertEqual(self.run_mock(plan,fake,clock=lambda:next(values)),1)
            fake.assert_not_called()

    def test_no_job_id_no_execution(self):
        with tempfile.TemporaryDirectory() as t:
            fake=Mock()
            with self.assertRaises(ValueError):worker.run_plan(self.temporary_plan(Path(t)),None,runner=fake)
            fake.assert_not_called()

    def test_hard_timeout_terminates_process_group(self):
        with tempfile.TemporaryDirectory() as t:
            proc=Mock(pid=123456789);proc.wait.side_effect=[subprocess.TimeoutExpired('MOCK',1),0]
            with patch.object(worker.subprocess,'Popen',return_value=proc),patch.object(worker.os,'killpg') as kill:
                result=worker.execute(['MOCK_ONLY'],1,Path(t)/'log')
            self.assertEqual(result['returncode'],124)
            kill.assert_called_once_with(proc.pid,signal.SIGTERM)
            self.assertIsNone(worker.ACTIVE)

    def test_hard_timeout_kill_after_grace(self):
        with tempfile.TemporaryDirectory() as t:
            proc=Mock(pid=123456789);proc.wait.side_effect=[subprocess.TimeoutExpired('MOCK',1),subprocess.TimeoutExpired('MOCK',120),0]
            with patch.object(worker.subprocess,'Popen',return_value=proc),patch.object(worker.os,'killpg') as kill:
                result=worker.execute(['MOCK_ONLY'],1,Path(t)/'log')
            self.assertEqual(result['returncode'],124)
            self.assertEqual([c.args[1] for c in kill.call_args_list],[signal.SIGTERM,signal.SIGKILL])

    def test_child_storage_error_stops(self):
        with tempfile.TemporaryDirectory() as t:
            def fake(argv,**options):
                options['stdout'].write(b'MOCK Disk quota exceeded\n');options['stdout'].flush()
                return Mock(pid=123456789,wait=Mock(return_value=1))
            with patch.object(worker.subprocess,'Popen',side_effect=fake):
                with self.assertRaises(OSError):worker.execute(['MOCK_ONLY'],1,Path(t)/'log')
    def test_main_arms_hard_shard_alarm(self):
        allocation={1:['G'+str(i) for i in range(7)]}
        with patch.object(sys,'argv',['worker','gp','01']),patch.dict(os.environ,{'SLURM_JOB_ID':'MOCK_JOB'}),patch.object(worker,'load_context',return_value=({}, {}, allocation)),patch.object(worker,'run_plan',return_value=0) as run,patch.object(worker.signal,'signal'),patch.object(worker.signal,'alarm') as alarm:
            self.assertEqual(worker.main(),0)
            run.assert_called_once()
            self.assertEqual([x.args[0] for x in alarm.call_args_list],[1209600,0])
    def test_bad_worker_pin_refused(self):
        config=json.loads(worker.CONFIG.read_text());config['worker_sha256']='wrong'
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'config.json';path.write_text(json.dumps(config))
            with patch.object(worker,'CONFIG',path):
                with self.assertRaises(ValueError):worker.load_context()

    def test_cli_plan_has_no_science_side_effect(self):
        before=sorted(str(p) for p in worker.SCIENCE.glob('gp_d5000_a16/*'))
        result=subprocess.run([sys.executable,'-B',worker.__file__,'gp','01','--plan'],capture_output=True,text=True,timeout=10)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(len(json.loads(result.stdout)['entries']),7)
        self.assertEqual(before,sorted(str(p) for p in worker.SCIENCE.glob('gp_d5000_a16/*')))

    def test_cli_login_execution_refused(self):
        env=dict(os.environ);env.pop('SLURM_JOB_ID',None)
        result=subprocess.run([sys.executable,'-B',worker.__file__,'gp','01'],env=env,capture_output=True,text=True,timeout=10)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('SLURM_JOB_ID absent',result.stderr)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt-root',type=Path,required=True)
    args=parser.parse_args()
    args.receipt_root.mkdir(parents=True,exist_ok=False)
    buffer=io.StringIO()
    result=unittest.TextTestRunner(stream=buffer,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(PackagingTests))
    print(buffer.getvalue())
    with (args.receipt_root/'unittest.log').open('x') as stream:stream.write(buffer.getvalue())
    worker.write_new(args.receipt_root/'test_receipt.json',{'tests':result.testsRun,'failures':len(result.failures),
        'errors':len(result.errors),'successful':result.wasSuccessful(),'worker_sha256':worker.sha(worker.__file__),
        'test_sha256':worker.sha(__file__),'config_sha256':worker.sha(worker.CONFIG),
        'production_manifest_sha256':worker.sha(worker.MANIFEST),'science_executions':0,'queue_submissions':0,
        'note':'actual scientific runner mocked; only subprocess plan/refusal and bash-n checks execute'})
    sys.exit(0 if result.wasSuccessful() else 1)

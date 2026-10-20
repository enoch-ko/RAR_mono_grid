"""Synthetic/mocked-only signed-rank tests: NO galaxy or force/probability solves."""
import copy
import io
import json
import math
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch
import warnings
import numpy as np
import signed_rank_policy_20261016 as policy
import run_signed_rank_20261016 as launcher


class SignedTests(unittest.TestCase):
    def baseline(self, nominal=(-1., 0.), nodes=((-1., 0.), (0., -1.)), weights=(.25, .75), **kwargs):
        return policy.SignedFrozenBaseline([1., 2.], [1., 1.], nominal, nodes, weights,
                                          integral=kwargs.pop('integral', lambda *a, **k: -.5), **kwargs)

    def evaluator(self, D=None):
        obj = policy.SignedProbabilityEvaluator.__new__(policy.SignedProbabilityEvaluator)
        obj.D = np.diag([1., -1.]) * 1e-12 if D is None else D
        obj.nuisance_terms = np.zeros((4, 2));obj.floor=1e-8
        obj.r=np.array([1., 2.]);obj.vobs=np.array([10., 20.]);obj.errv=np.ones(2)
        obj.vgas=np.zeros(2);obj.vbul=np.zeros(2)
        obj.gas_nodes=np.array([.8,1.2]);obj.gas_weights=np.array([.25,.75])
        obj.bulge_nodes=np.array([.5,.9]);obj.bulge_weights=np.array([.4,.6])
        obj.args=SimpleNamespace(gas_ml=1.,bulge_ml=.7,marginalise_gas_bulge=True,
                                 prob_grid_size=400,zpadding=8.,tie_atol=0.,max_tie_permutations=720)
        obj._ordered_logp=Mock(return_value=-1.)
        return obj

    def test_parent_pins(self):
        policy.verify_parents()

    def test_syntax_all_new(self):
        for name in (*launcher.NEW_SOURCES,Path(__file__).name):
            compile((policy.ROOT/name).read_text(),name,'exec')

    def test_signed_nominal_and_nodes(self):
        b=self.baseline();self.assertAlmostEqual(b.score([1.,2.]),-.5)
        np.testing.assert_array_equal(b.nodes,[[-1,0],[0,-1]])
        self.assertFalse(b.nodes.flags.writeable)

    def test_zero_ties_average_not_zero_probability(self):
        fake=lambda v,*a,**k: math.log(.2 if v[0]==1 else .8)
        b=self.baseline(nominal=(0.,0.),nodes=((0.,0.),),weights=(1.,),integral=fake)
        self.assertEqual(len(b.orders[0]),2)
        self.assertAlmostEqual(math.exp(b.score([1.,2.])),.5)

    def test_every_node_weight_retained(self):
        fake=lambda v,*a,**k: math.log(.2 if v[0]==2 else .8)
        b=self.baseline(integral=fake)
        np.testing.assert_array_equal(b.weights,[.25,.75])
        self.assertAlmostEqual(math.exp(b.score([1.,2.])),.25*.2+.75*.8)

    def test_positive_reference_matches_parent_mocked_integral(self):
        args=([1.,2.],[1.,1.],[1.,2.],[[1.,2.],[2.,1.]],[.25,.75])
        fake=lambda v,*a,**k: math.log(.2 if v[0]==2 else .8)
        old=policy.common.FrozenBaseline(*args,integral=fake)
        new=policy.SignedFrozenBaseline(*args,integral=fake)
        self.assertEqual(old.orders,new.orders)
        self.assertEqual(old.score([1.,2.]),new.score([1.,2.]))

    def test_nonfinite_reference_refused(self):
        for value in (np.nan,np.inf,-np.inf):
            with self.subTest(value=value),self.assertRaises(ValueError):
                self.baseline(nodes=((value,0.),(0.,-1.)))

    def test_invalid_weights_refused(self):
        for weights in ((0.,1.),(-1.,2.),(.2,.2),(np.nan,1.)):
            with self.subTest(weights=weights),self.assertRaises(ValueError):
                self.baseline(weights=weights)

    def test_reference_positive_radius_error_still_required(self):
        for r,e in [([0.,2.],[1.,1.]),([1.,2.],[0.,1.])]:
            with self.assertRaises(ValueError):
                policy.SignedFrozenBaseline(r,e,[-1.,0.],[[-1.,0.]],[1.],integral=Mock())

    def test_tie_cap_still_refused(self):
        with self.assertRaises(ValueError):
            self.baseline(nominal=(0.,0.),nodes=((0.,0.),),weights=(1.,),
                          settings=policy.common.Settings(max_tie_permutations=1))

    def test_nonfinite_probability_still_refused(self):
        for logp in (np.nan,np.inf,-np.inf,.1):
            with self.subTest(logp=logp),self.assertRaises(ValueError):
                self.baseline(integral=lambda *a,**k:logp).score([1.,2.])

    def test_probability_warning_still_fatal(self):
        def fake(*a,**k):warnings.warn('MOCK numerical warning');return -.5
        with self.assertRaises(UserWarning):self.baseline(integral=fake).score([1.,2.])

    def test_isotonic_accepts_signed_abscissa(self):
        fit=policy.common.isotonic_reference([-1.,0.,0.,2.],[1.,2.,3.,4.],[4.,3.,2.,1.],[1.]*4)
        theta=np.array(fit['theta']);order=fit['ascending_order']
        self.assertTrue(np.all(theta>=0));self.assertTrue(np.all(np.diff(theta[order])>=0))
        self.assertEqual(theta[1],theta[2])

    def test_disc_interval_no_gbar_constraint(self):
        mean=np.array([.5,.5]);direction=np.array([1.,-1.]);D=-np.eye(2)
        old=policy.gp.physical_interval(mean,direction,.01,10.,D,np.zeros((1,2)),1e-8)
        self.assertIsNone(old)
        self.assertEqual(policy.physical_interval(mean,direction,.01,10.,D,np.zeros((1,2)),1e-8),(0.,.49))

    def test_interval_positive_disc_bounds(self):
        for lower in (0.,-.1,np.nan):
            with self.assertRaises(ValueError):policy.physical_interval([.5],[1.],lower,10.)
        self.assertIsNone(policy.physical_interval([np.nan],[1.],.01,10.))

    def test_signed_gp_evaluation_uses_all_nodes(self):
        obj=self.evaluator();result=obj.evaluate([.5,.5])
        self.assertTrue(result['valid']);self.assertEqual(result['nonpositive_nuisance_entries'],4)
        self.assertFalse(result['positive_net_acceleration_at_all_nodes'])
        self.assertEqual(obj._ordered_logp.call_count,5)
        self.assertTrue(all(np.any(c.args[0]<0) for c in obj._ordered_logp.call_args_list))
        self.assertAlmostEqual(result['marginal_logp'],-1.)

    def test_zero_gp_values_not_clipped(self):
        obj=self.evaluator(np.zeros((2,2)));result=obj.evaluate([.5,.5])
        self.assertTrue(result['valid']);self.assertEqual(result['nonpositive_nuisance_entries'],8)
        self.assertTrue(all(np.all(c.args[0]==0) for c in obj._ordered_logp.call_args_list))

    def test_gp_disc_and_nonfinite_force_safeguards(self):
        for xi in ([0.,.5],[-1.,.5],[np.nan,.5]):
            obj=self.evaluator();self.assertFalse(obj.evaluate(xi)['valid']);obj._ordered_logp.assert_not_called()
        obj=self.evaluator(np.full((2,2),np.nan))
        self.assertFalse(obj.evaluate([.5,.5])['valid']);obj._ordered_logp.assert_not_called()

    def test_gp_numerical_failure_and_warning_refused(self):
        obj=self.evaluator()
        for result in [((np.nan,False,'MOCK',[]),(-1.,True,'',[])),
                       ((-1.,True,'',['MOCK warning']),(-1.,True,'',[]))]:
            with patch.object(policy.gp.mono,'evaluate_disk_term',return_value=result):
                self.assertFalse(obj.evaluate([.5,.5])['valid'])

    def test_display_mask_does_not_change_scoring_array(self):
        a=np.array([-1.,0.,2.,np.nan]);before=a.copy()
        np.testing.assert_array_equal(policy.log_axis_mask(a),[False,False,True,False])
        np.testing.assert_array_equal(a,before)

    def test_context_is_versioned_not_mutated(self):
        old={'invariants':{'parser':'MOCK'}};before=copy.deepcopy(old)
        new=policy.signed_context(old)
        self.assertEqual(old,before);self.assertNotEqual(new['invariants'],old['invariants'])
        self.assertEqual(new['invariants']['signed_rank_policy']['id'],policy.POLICY_ID)

    def test_active_manifest_refused(self):
        active=json.loads((policy.ROOT/'outputs/science_20261016/manifest.json').read_text())
        with self.assertRaises(ValueError):launcher.validate_manifest(active)

    def test_new_manifest_requires_namespace_and_new_pins(self):
        d={'signed_rank_policy':policy.POLICY_ID,'pilot_output_root':str(launcher.NEW_ROOT),
           'source_sha256':{n:policy.sha(policy.ROOT/n) for n in launcher.NEW_SOURCES}}
        launcher.validate_manifest(d)
        wrong=copy.deepcopy(d);wrong['pilot_output_root']=str(policy.ROOT/'outputs/science_20261016')
        with self.assertRaises(ValueError):launcher.validate_manifest(wrong)
        wrong=copy.deepcopy(d);wrong['source_sha256']={}
        with self.assertRaises(ValueError):launcher.validate_manifest(wrong)

    def test_describe_has_no_science(self):
        with patch('sys.stdout',new_callable=io.StringIO) as out:launcher.main(['--describe'])
        d=json.loads(out.getvalue());self.assertEqual(d['science_evaluations'],0)
        self.assertEqual(d['queue_submissions'],0)

    def test_bad_parent_pin_refused(self):
        with patch.dict(policy.PARENT_PINS,{'common_isotonic_reference.py':'WRONG'}):
            with self.assertRaises(ValueError):policy.verify_parents()


    def test_reference_install_context_and_save_are_versioned_mocked(self):
        core=policy.common
        fake_run=Mock(side_effect=lambda velocity,scorer,**kw: {'context':kw['context'],'status':'MOCK'})
        fake_save=Mock(return_value='MOCK_SAVED')
        with patch.object(core,'_signed_rank_policy_installed',False,create=True),patch.multiple(
                core,FrozenBaseline=core.FrozenBaseline,REFERENCE_ID=core.REFERENCE_ID,
                run_reference=fake_run,save_reference=fake_save):
            policy.install_reference()
            self.assertIs(core.FrozenBaseline,policy.SignedFrozenBaseline)
            context={'invariants':{'MOCK':True}}
            result=core.run_reference([1.,2.],self.baseline(),B=8,seed=1,context=context)
            self.assertNotIn('signed_rank_policy',context)
            self.assertEqual(result['context']['invariants']['signed_rank_policy']['id'],policy.POLICY_ID)
            self.assertEqual(result['nonpositive_nuisance_entries'],4)
            self.assertEqual(core.save_reference(result,'MOCK_PATH'),'MOCK_SAVED')
            self.assertIn('signed_rank_policy_20261016.py',result['driver_source_sha256'])
            fake_save.assert_called_once()

    def test_gp_install_changes_only_mocked_process_bindings(self):
        g=policy.gp;fake=Mock(return_value={'settings':{},'physical_domain':{}})
        with patch.object(g,'_signed_rank_policy_installed',False,create=True),patch.multiple(
                g,physical_interval=g.physical_interval,ProbabilityEvaluator=g.ProbabilityEvaluator,
                SOURCE_FILES=g.SOURCE_FILES,galaxy_run=fake):
            policy.install_gp();args=SimpleNamespace()
            result=g.galaxy_run('SYNTHETIC',args)
            self.assertEqual(args.signed_rank_policy,policy.POLICY_ID)
            self.assertIsNone(result['physical_domain']['floor_m_s2'])
            self.assertIsNone(result['settings']['signed_gbar_floor_relative'])
            self.assertIn('signed_rank_policy_20261016.py',g.SOURCE_FILES)
            self.assertIs(g.ProbabilityEvaluator,policy.SignedProbabilityEvaluator)
            with self.assertRaises(ValueError):g.galaxy_run('SYNTHETIC',args,'OLD_QP.json')
            fake.assert_called_once()


if __name__=='__main__':
    unittest.main(verbosity=2)

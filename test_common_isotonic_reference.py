"""Standalone bounded synthetic tests: python3 -B -m unittest -v test_common_isotonic_reference."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
import warnings
import numpy as np
from scipy.optimize import minimize
import common_isotonic_reference as c


class CommonReferenceTests(unittest.TestCase):
    def scorer(self, integral=None):
        return c.FrozenBaseline([1., 2.], [1., 1.5], [2., 1.],
                                [[2., 1.], [2., 1.]], [.3, .7], integral=integral)

    def context(self):
        return dict(galaxy="synthetic", parser_id=c.PARSER_ID,
                    input_sha256={"fixture": "synthetic"}, normalization_id="fixture_only",
                    nuisance_model="two_duplicate_nodes_fixture")

    def test_header_gas_disk_sentinel(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/"fixture.dat"
            path.write_text("# Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul\n"
                            "0 5 1 7 70 0 1 0\n1 8 1 -3 31 2 1 0\n2 9 1 4 42 3 1 0\n")
            data, meta=c.read_rotmod(path, (0, 2))
            np.testing.assert_array_equal(data["Vgas"], [-3, 4])
            np.testing.assert_array_equal(data["Vdisk"], [31, 42])
            self.assertEqual(meta["mask"], [False, True, True])
            self.assertEqual(meta["parser_id"], c.PARSER_ID)
            path.write_text(path.read_text().replace("Vgas Vdisk", "Vdisk Vgas"))
            with self.assertRaises(ValueError): c.read_rotmod(path, (0, 2))

    def test_isotonic_ties_restore_nonnegative(self):
        fit=c.isotonic_reference([2, 1, 1, 3], [1]*4, [4, 3, 1, -2], [1]*4)
        np.testing.assert_allclose(fit["theta"], [1.5]*4)
        self.assertEqual(fit["exact_tie_groups"], [[1, 2]])
        fit=c.isotonic_reference([1, 2], [1, 4], [-2, -1], [1, 2])
        np.testing.assert_array_equal(fit["theta"], [0, 0])

    def test_pava_independent_constrained_minimum(self):
        rng=np.random.default_rng(77)
        for _ in range(12):
            x=np.arange(6.); r=np.arange(1.,7.); e=rng.uniform(.5,2.,6); v=rng.normal(1,3,6)
            fit=np.array(c.isotonic_reference(x,r,v,e)["theta"])
            y=v/np.sqrt(r); w=r/e**2
            result=minimize(lambda t: np.sum(w*(t-y)**2), np.maximum.accumulate(np.maximum(y,0)),
                            jac=lambda t:2*w*(t-y), method="SLSQP", bounds=[(0,None)]*6,
                            constraints={"type":"ineq","fun":lambda t:np.diff(t)},
                            options={"ftol":1e-11,"maxiter":500})
            self.assertTrue(result.success, result.message)
            np.testing.assert_allclose(fit,result.x,atol=2e-6,rtol=2e-6)

    def test_per_draw_cache_invalidation(self):
        calls=[]
        def integral(v,e,r,**kw):
            calls.append(v.copy()); return -float(np.sum(v*v))
        scorer=self.scorer(integral)
        self.assertEqual(scorer.score([1,2]),-5)
        self.assertEqual(scorer.score([2,3]),-13)
        self.assertEqual(len(calls),2)

    def test_warning_nonfinite_and_floor_rejected(self):
        def warning(*a,**k): warnings.warn("fixture warning",RuntimeWarning); return -1
        with self.assertRaises(RuntimeWarning): self.scorer(warning).score([1,2])
        with self.assertRaises(ValueError): self.scorer(lambda *a,**k:-np.inf).score([1,2])
        with self.assertRaises(ValueError):
            c.FrozenBaseline([1,2],[1,1],[2,1],[[2,1e-20]],[1])
        with self.assertRaises(ValueError):
            c.FrozenBaseline(np.ones(7),np.ones(7),np.ones(7),[np.ones(7)],[1])

    def test_real_backend_tied_and_distinct_equivalence(self):
        from probability_recursion import boundary_log_probability_for_gbar
        from utils import log_weighted_average
        r=np.array([1.,2.,3.]);e=np.array([.5,.8,1.]);v=np.array([2.,-1.,3.])
        nodes=np.array([[3.,2.,1.],[2.,2.,1.]]);w=np.array([.4,.6])
        frozen=c.FrozenBaseline(r,e,nodes[0],nodes,w)
        direct=[boundary_log_probability_for_gbar(g,v,e,r,ngrid=400,zpadding=8.,
                    tie_atol=0.,max_tie_permutations=720,log_tolerance=2e-5,max_refinements=8)
                for g in nodes]
        self.assertAlmostEqual(frozen.score(v),log_weighted_average(np.array(direct),w),places=12)

    def test_mock_replay_signed_velocities_and_context(self):
        scorer=self.scorer(lambda v,e,r,**kw:-float(np.sum(v*v)))
        a=c.run_reference([-5,-5],scorer,B=20,seed=11,context=self.context(),save_vectors=True)
        b=c.run_reference([-5,-5],scorer,B=20,seed=11,context=self.context(),save_vectors=True)
        self.assertEqual(a,b); self.assertEqual(a["status"],"complete")
        np.testing.assert_array_equal(a["mock_velocity"],
            np.random.Generator(np.random.PCG64(11)).standard_normal((20,2))*scorer.e)
        self.assertGreater(a["negative_velocity_draws"],0)
        self.assertEqual(a["quoted_error"],[1.,1.5])
        json.dumps(a,allow_nan=False)
        self.assertIn("q_common_ref",c.annotate(a,-1,a["context"]))
        bad=copy.deepcopy(a["context"]);bad["parser_id"]="legacy_gas4"
        with self.assertRaises(ValueError): c.annotate(a,-1,bad)
        bad=copy.deepcopy(a);bad["context"]["galaxy"]="other"
        with self.assertRaises(ValueError): c.annotate(bad,-1,bad["context"])

    def test_failures_do_not_shrink_denominator(self):
        scorer=self.scorer(lambda *a,**k:float("nan"))
        out=c.run_reference([1,2],scorer,B=3,seed=1,context=self.context())
        self.assertEqual(out["status"],"unusable");self.assertEqual(out["B"],3)
        self.assertEqual(out["mock_logp"],[None]*3);self.assertEqual(len(out["failures"]),4)
        with self.assertRaises(ValueError): c.annotate(out,-1,out["context"])

    def test_inclusive_cdf_cutoffs_and_uncertainty(self):
        for values in [[-4,-2,-2,-1],list(-np.arange(1,101,dtype=float))]:
            cdf=c.EmpiricalCDF(values)
            for alpha in [0,.01,.05,.07,.14,.28,.29,.5,1]:
                cut=cdf.cutoff(alpha)
                for t in [-101,-4,-3,-2,-1,0]:
                    expected=cdf.evaluate(t)["q_common_ref"] >= alpha
                    got=cut["accept_all_valid_scores"] or t>=cut["logp_cutoff"]
                    self.assertEqual(expected,got,(alpha,t))
            for t in [-101,-2,0]:
                result=cdf.evaluate(t);lo,hi=result["conditional_mc_wilson95"]
                self.assertLessEqual(lo,result["q_common_ref"]+1e-15)
                self.assertGreaterEqual(hi,result["q_common_ref"]-1e-15)
        self.assertEqual(c.EmpiricalCDF([-4,-2,-2,-1]).evaluate(-2)["q_common_ref"],.75)
        for bad in [[-np.inf],[np.nan],[],[.1]]:
            with self.assertRaises(ValueError): c.EmpiricalCDF(bad)

    def test_fixed_baseline_uses_signed_gas_and_casertano_input(self):
        from utils import vel2acc
        data={"Rad":np.array([1.,2.]),"errV":np.array([1.,1.]),
              "Vgas":np.array([-2.,3.]),"Vbul":np.array([1.,2.]),
              "Vdisk":np.array([999.,999.])}
        w=np.array([100.,100.])
        scorer=c.prepare_sparc_baseline(data,w)
        np.testing.assert_allclose(scorer.nominal,vel2acc(.5*w+[-4,9]+.7*data["Vbul"]**2,data["Rad"]))
        self.assertEqual(scorer.nodes.shape,(81,2))
        self.assertAlmostEqual(scorer.weights.sum(),1.)

    def test_bounded_exclusive_artifact_save(self):
        out=c.run_reference([1,2],self.scorer(lambda *a,**k:-1.),B=2,seed=1,context=self.context())
        self.assertTrue(out["degenerate_point_mass"])
        self.assertIn("common_isotonic_reference.py",out["context"]["source_sha256"])
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"reference.json"
            with self.assertRaises(ValueError): c.save_reference(out,p,max_bytes=1)
            self.assertFalse(p.exists())
            c.save_reference(out,p)
            self.assertEqual(json.loads(p.read_text()),out)
            with self.assertRaises(FileExistsError): c.save_reference(out,p)

    def test_full_GH_adapter_matches_correct_grid_marginalizer(self):
        from functools import partial
        from types import SimpleNamespace
        from unittest.mock import patch
        import mono_SPARC_RAR as mono
        from probability_recursion import boundary_log_probability_for_gbar
        from utils import lognormal_quadrature
        data={"Rad":np.array([1.,2.,3.]),"errV":np.array([.5,.8,1.]),
              "Vobs":np.array([2.,-1.,3.]),"Vgas":np.array([-2.,3.,1.]),
              "Vbul":np.array([1.,2.,1.])}
        w=np.array([100.,80.,60.]); frozen=c.prepare_sparc_baseline(data,w)
        gn,gw=lognormal_quadrature(1.,.04,9);bn,bw=lognormal_quadrature(.7,.1,9)
        np.testing.assert_allclose(frozen.weights,np.outer(gw,bw).ravel(),rtol=1e-14)
        args=SimpleNamespace(gas_ml=1.,bulge_ml=.7,marginalise_gas_bulge=True,
            prob_grid_size=400,zpadding=8.,tie_atol=0.,max_tie_permutations=720)
        with patch.object(mono,"log_monotonic_probability_for_gbar",partial(
                boundary_log_probability_for_gbar,log_tolerance=2e-5,max_refinements=8)):
            cond,marg=mono.evaluate_disk_term(.5*w,data["Vgas"],data["Vbul"],data["Rad"],
                data["Vobs"],data["errV"],args,gn,gw,bn,bw)
        self.assertTrue(marg[1],marg[2]);self.assertFalse(marg[3])
        self.assertAlmostEqual(frozen.score(data["Vobs"]),marg[0],places=12)

    def test_different_candidate_profile_same_independent_invariants(self):
        from dataclasses import asdict
        import common_reference_driver as driver
        scorer=self.scorer(lambda *a,**k:-1.)
        def independently_construct_evidence():
            return {"galaxy":"synthetic","parser_id":c.PARSER_ID,
                "input_sha256":{"raw_fixture":"sha-fixture"},"normalization_id":"synthetic",
                "retained_mask":[True,False,True],"nuisance_model":"synthetic",
                "kinematics_sha256":c.digest({"r":[1.,2.],"v":[1.,2.],"e":[1.,1.5]}),
                "probability_settings":asdict(scorer.settings),
                "force_model":{"geometry":{"zd_kpc":.196},"interpolation":"signed_force",
                               "force_source_sha256":{"fixture":"force-v1"}}}
        context=self.context();context["invariants"]=independently_construct_evidence()
        reference=c.run_reference([1.,2.],scorer,B=3,seed=2,context=context)
        doc={"schema":"common-reference-candidates-v1", "invariants":independently_construct_evidence(),
             "deviation_metric":{"id":"synthetic_fixed_J_fixture"},
             "parser_source_sha256":reference["context"]["source_sha256"]["common_isotonic_reference.py"],
             "records":[{"valid":True,"probability_kind":"marginalised_gas_bulge",
                         "disc_ml":.8,"disk_slope":.03,"log_probability":-.5,"J_ref":4.}]}
        output=driver.annotate_document(reference,doc)
        self.assertEqual(output["records"][0]["disc_ml"],.8)
        self.assertEqual(output["records"][0]["J_ref"],4.)
        self.assertEqual(output["records"][0]["q_common_ref"],1.)
        self.assertEqual(output["minimum_found_J_record_index"],0)
        bad=copy.deepcopy(doc);bad["invariants"]["force_model"]["geometry"]["zd_kpc"]*=2
        with self.assertRaises(ValueError):driver.annotate_document(reference,bad)
        bad=copy.deepcopy(doc);bad["invariants"]["force_model"]["interpolation"]="velocity"
        with self.assertRaises(ValueError):driver.annotate_document(reference,bad)
        bad=copy.deepcopy(doc);bad.pop("deviation_metric")
        with self.assertRaises(ValueError):driver.annotate_document(reference,bad)
        bad=copy.deepcopy(doc);bad["invariants"]["galaxy"]="wrong-name"
        with self.assertRaises(ValueError):driver.annotate_document(reference,bad)
        bad=copy.deepcopy(doc);bad["invariants"]["retained_mask"]=[True,True,False]
        with self.assertRaises(ValueError):driver.annotate_document(reference,bad)
        bad=copy.deepcopy(doc);bad["provenance"]={"source_sha256":{"inverse_rar_diagnostic.py":driver.BAD_INVERSE_HASH}}
        with self.assertRaisesRegex(ValueError,"wrong-gas-column"):driver.annotate_document(reference,bad)
        with self.assertRaises(ValueError):driver.annotate_document(reference,{"records":[]})

    def test_timeout_is_not_hidden_as_mock_failure(self):
        def integral(*a,**kw): raise TimeoutError("fixture cap")
        with self.assertRaises(TimeoutError):
            c.run_reference([1,2],self.scorer(integral),B=2,seed=1,context=self.context())

    def test_units_and_radial_index_sentinel(self):
        from utils import vel2acc,compute_A_B
        self.assertAlmostEqual(float(vel2acc(1.,1.)),1e6/3.086e19,places=25)
        luminosity=compute_A_B(np.array([0.,1.]),np.array([2.,2.])*1e6)[0]
        self.assertAlmostEqual(luminosity/1e6,2*np.pi,places=12)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/"rot.dat"
            p.write_text("# Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul\n"
                ".2 11 1 -2 22 3 4 5\n1 33 2 -6 66 7 8 9\n3 44 3 -10 100 11 12 13\n")
            values,meta=c.read_rotmod(p,(.5,2.))
            for key,value in {"Rad":1,"Vobs":33,"errV":2,"Vgas":-6,"Vdisk":66,"Vbul":7}.items():
                self.assertEqual(values[key].tolist(),[value])
            self.assertEqual(meta["mask"],[False,True,False])

    def test_driver_fixture_no_force_dry_run_and_one_solve_build(self):
        import common_reference_driver as driver
        import pandas as pd
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);galaxy="Fixture"
            (root/(galaxy+".dens")).write_text("# Rad[kpc] SBdisk[Lsun/pc^2] SBbulge[Lsun/pc^2]\n"
                "0\t3\t1\n1\t2\t.5\n2\t1\t.2\n")
            (root/(galaxy+"_rotmod.dat")).write_text("# Distance = 8 Mpc\n"
                "# Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul\n"
                "0 7 1 1 99 1 1 1\n1 8 1 -3 31 2 1 1\n2 9 1 4 42 3 1 1\n")
            (root/"table.txt").write_text("synthetic table fixture only")
            norm={"policy":"disc_to_catalogue_total_minus_published_bulge_v1","galaxies":{galaxy:{
                "Ltotal_1e9":2.,"Lbulge_1e9":.5,"Ldisk_1e9":1.5,"kappa":.75,"reference_distance_Mpc":8.}}}
            (root/"norm.json").write_text(json.dumps(norm))
            args=driver.parser().parse_args(["build","--galaxy",galaxy,"--normalization",str(root/"norm.json"),
                "--B","2","--seed","3","--output",str(root/"result.json"),"--data-root",str(root),
                "--dens-root",str(root),"--table",str(root/"table.txt"),"--dry-run"])
            table=pd.DataFrame([{"Galaxy":galaxy,"D":8.,"L":2.,"Rdisk":1.,"Inc":60.}])
            calls=[]
            def synthetic_solve(rad,sb,zd,r):
                calls.append((rad.copy(),r.copy()))
                return np.full(len(r),100.),[],""
            with patch.object(driver.utils,"load_sparc_table",return_value=table):
                plan=driver.build(args,solve=synthetic_solve)
                self.assertEqual(plan["status"],"prepared_no_force_no_draws")
                self.assertEqual(calls,[]);self.assertFalse(args.output.exists())
                args.dry_run=False
                result=driver.build(args,solve=synthetic_solve)
                self.assertEqual(result["status"],"complete")
                self.assertEqual(len(calls),1)
                saved=json.loads(args.output.read_text())
                self.assertEqual(saved["B"],2)
                self.assertEqual(saved["context"]["invariants"]["retained_mask"],[False,True,True])
                with self.assertRaises(FileExistsError):driver.build(args,solve=synthetic_solve)
                bad=dict(norm);bad["galaxies"]={"Other":norm["galaxies"][galaxy]}
                (root/"norm.json").write_text(json.dumps(bad))
                with self.assertRaises(KeyError):driver.load_inputs(galaxy,root/"norm.json",root,root,root/"table.txt")


if __name__ == "__main__":
    unittest.main()

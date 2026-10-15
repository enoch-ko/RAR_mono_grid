"""Independent bounded checks; no production galaxy search or fit."""
import copy
import math
import unittest
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import norm
import report_mg_compactness as reporter
import mock_fixed_profile_N as mocks


def pair_probability(v, e, r):
    """Independent full folded-normal integral P(|V1|/sqrt(r1)>|V2|/sqrt(r2))."""
    mu = np.asarray(v)/np.sqrt(r)
    sd = np.asarray(e)/np.sqrt(r)
    def integrand(x):
        density = (norm.pdf((x-mu[0])/sd[0])+norm.pdf((x+mu[0])/sd[0]))/sd[0]
        cdf = ndtr((x-mu[1])/sd[1])-ndtr((-x-mu[1])/sd[1])
        return density*cdf
    value, error = quad(integrand,0.,np.inf,epsabs=1e-11,epsrel=1e-10,limit=300)
    if error > 1e-8: raise AssertionError("independent pair quadrature inaccurate")
    return value


class StandaloneDiagnosticTests(unittest.TestCase):
    def test_mass_availability_support_and_retained_adapter(self):
        reporter.self_test()

    def test_nested_null_and_probability_wrapper(self):
        mocks.self_test()

    def test_one_bin(self):
        lp, _, _ = mocks.fixed_probability([1.],[0.],[0.],[-.5],[1.],[2.],ngrid=64)
        self.assertEqual(lp,0.)

    def test_independent_unequal_pair(self):
        self.check_pair([20.,12.],[3.,2.],[2.,1.])

    def test_independent_low_snr_pair_both_signs(self):
        self.check_pair([.3,-.8],[1.2,.9],[1.,3.])

    def check_pair(self,v,e,r):
        arrays=[np.asarray(x) for x in (v,e,r)]
        p=pair_probability(*arrays)
        swapped=pair_probability(*(a[::-1] for a in arrays))
        self.assertAlmostEqual(p+swapped,1.,places=9)
        lp,_,_=mocks.fixed_probability([2.,1.],[0.,0.],[0.,0.],v,e,r,ngrid=101)
        sp,_,_=mocks.fixed_probability([2.,1.],[0.,0.],[0.,0.],v[::-1],e[::-1],r[::-1],ngrid=101)
        self.assertLess(abs(math.exp(lp)-p),2e-5)
        self.assertLess(abs(math.exp(sp)-swapped),2e-5)
        print("PAIR_CHECK",v,"independent",p,"recursion",math.exp(lp),"swap_sum",math.exp(lp)+math.exp(sp))

    def test_iid_three_accelerations(self):
        lp,_,_=mocks.fixed_probability([3.,2.,1.],[0.,0.,0.],[0.,0.,0.],
                                       [5.,10.,15.],[1.,2.,3.],[1.,4.,9.],ngrid=101)
        self.assertLess(abs(math.exp(lp)-1/6),2e-5)
        print("IID_N3",math.exp(lp),"expected",1/6)

    def test_model_not_defined_preserves_mixed_proxy(self):
        row=self.row()
        reporter.apply_model_specs([row],None,1.,1.2e-10)
        self.assertEqual(row["model_mass_status"],"not_defined")
        self.assertIsNone(row["eta_model"])
        self.assertEqual(row["eta_mixed_proxy"],.9)
        self.assertIsNone(row["eta_total"])

    def test_complete_declared_model(self):
        row=self.row();spec=self.spec()
        reporter.apply_model_specs([row],spec,1.,1.2e-10)
        self.assertAlmostEqual(row["Mbar_model_1e9"],.6)
        self.assertGreater(row["eta_model"],0.)
        self.assertEqual(row["force_model_match"],"unverified_supplied_model")
        self.assertIsNone(row["eta_total"])

    def test_missing_model_bulge_is_not_zero(self):
        row=self.row();spec=self.spec();del spec["entries"][0]["components"]["bulge"]
        reporter.apply_model_specs([row],spec,1.,1.2e-10)
        self.assertEqual(row["component_completeness"],"missing:bulge")
        self.assertIsNone(row["Mbar_model_1e9"])

    def test_model_boundaries_and_disc_conflict_rejected(self):
        spec=self.spec();del spec["entries"][0]["components"]["gas"]["outer_boundary"]
        with self.assertRaises(ValueError): reporter.apply_model_specs([self.row()],spec,1.,1.2e-10)
        spec=self.spec();spec["entries"][0]["components"]["disk"]["mass_1e9"] = .9
        with self.assertRaises(ValueError): reporter.apply_model_specs([self.row()],spec,1.,1.2e-10)

    @staticmethod
    def row():
        return dict(profile_id="synthetic",convention="solver_disk_normalization",support="original_photometric",
                    M_disk_support_1e9=.5,eta_total=None,eta_mixed_proxy=.9)

    @classmethod
    def spec(cls):
        row=cls.row();entry={k:row[k] for k in ("profile_id","convention","support")}
        entry.update(full_model_defined=True,model_convention="synthetic_test_only",
                     components={name:dict(mass_1e9=mass,units="1e9_Msun",normalization="explicit synthetic",
                                           source="synthetic test, not a galaxy",inner_boundary="regular centre",
                                           outer_boundary="zero outside declared synthetic support")
                                 for name,mass in (("disk",.5),("gas",.1),("bulge",0.))})
        return dict(entries=[entry])


if __name__ == "__main__":
    unittest.main()

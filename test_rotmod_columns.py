"""Regression tests for isolated rotmod-column corrections; no real-galaxy fits."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch
import numpy as np
import pandas as pd
import sparc_rotmod_columns as columns
import inverse_rar_diagnostic_corrected as inverse
import plot_inverse_rar_comparison_corrected as plotter
import gp_rar_search_retained_corrected as gp


class RotmodColumnTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.rot=self.root/'Fixture_rotmod.dat'
        self.rot.write_text('# Distance = 8 Mpc\n# Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul\n'
            '.05 5 1 3 30 1 11 111\n1 8 2 -4 40 2 22 222\n2 9 3 5 50 3 33 333\n')
        (self.root/'Fixture.dens').write_text('# Rad SBdisk SBbulge\n.1 3 0\n1 2 0\n2 1 0\n')
        self.table=pd.DataFrame([{'Galaxy':'Fixture','Rdisk':1.}])
        self.rad=np.array([.1,1.,2.]);self.sb=np.array([3.,2.,1.])

    def tearDown(self):self.tmp.cleanup()

    def load(self):
        with patch.object(inverse,'load_sparc_table',return_value=self.table),patch.object(
                inverse,'get_SBdisk',return_value=(self.rad,self.sb)):
            return inverse.galaxy_data('Fixture',self.root,self.root)

    def test_named_columns_and_joint_radius_mask(self):
        rad,sb,r,v,e,gas,bulge,zd,excluded=self.load()
        np.testing.assert_array_equal(r,[1,2]);np.testing.assert_array_equal(v,[8,9])
        np.testing.assert_array_equal(e,[2,3]);np.testing.assert_array_equal(gas,[-4,5])
        np.testing.assert_array_equal(bulge,[2,3]);self.assertEqual(excluded,1)
        np.testing.assert_array_equal(rad,self.rad);self.assertAlmostEqual(zd,.196)

    def test_disc_reference_uses_disc_not_gas(self):
        raw=columns.read_rotmod(self.rot);mask=raw[:,0]>=.1
        disc=columns.disc_reference_velocity(raw,mask,.5)
        np.testing.assert_allclose(disc,np.array([40,50])*np.sqrt(.5))
        self.assertFalse(np.allclose(disc,np.array([-4,5])*np.sqrt(.5)))

    def test_reject_swapped_or_missing_header(self):
        original=self.rot.read_text()
        for text in [original.replace('Vgas Vdisk','Vdisk Vgas'),original.split('\n',2)[2]]:
            self.rot.write_text(text)
            with self.assertRaises(ValueError):columns.read_rotmod(self.rot)

    def test_provenance_has_semantic_indices_and_source_hash(self):
        proof=columns.parser_provenance()
        self.assertEqual(proof['columns']['Vgas'],3);self.assertEqual(proof['columns']['Vdisk'],4)
        self.assertEqual(proof['source_sha256'],hashlib.sha256(Path(columns.__file__).read_bytes()).hexdigest())
        self.assertIs(gp.inverse,inverse)
        self.assertIn('inverse_rar_diagnostic_corrected.py',gp.SOURCE_FILES)
        self.assertIn('sparc_rotmod_columns.py',gp.SOURCE_FILES)

    def test_QP_constraints_receive_signed_true_gas(self):
        values=self.load();seen=[]
        def fake_fit(D,c,observed,knots,args):
            seen.append(c.copy());return {'status':'infeasible'},[],10.
        args=SimpleNamespace(galaxy='Fixture',data_root=self.root,dens_root=self.root,
            validation_root=self.root,knots_from_dens=True,length_fraction=.3,
            min_ml=.01,max_ml=10.,gas_ml=1.,bulge_ml=.7)
        with patch.object(inverse,'galaxy_data',return_value=values),patch.object(inverse,'cache_gate',return_value={}),             patch.object(inverse,'build_response',return_value=(np.eye(3),np.zeros((2,3)),[])),             patch.object(inverse,'direct_check',return_value={}),patch.object(inverse,'solve_with_retry',side_effect=fake_fit):
            result=inverse.analyze(args)
        expected=inverse.vel2acc(np.array([-16,25])+.7*np.array([4,9]),np.array([1,2]))
        np.testing.assert_allclose(seen[0],expected)
        self.assertEqual(result['provenance']['input_parser'],columns.parser_provenance())

    def test_reject_legacy_QP_reference(self):
        path=self.root/'legacy.json';path.write_text(json.dumps({'provenance':{'script_sha256':'old'}}))
        with self.assertRaisesRegex(ValueError,'corrected gas-column'):
            gp.reference_profile(path,'Fixture',self.rad,.01,10.,None,{})

    def test_corrected_QP_reference_can_be_imported(self):
        args=SimpleNamespace(length_fraction=.3,gas_ml=1.,bulge_ml=.7,gbar_floor_relative=1e-8)
        data={'galaxy':'Fixture','input_sha256':{'fixture':'hash'},'knots_kpc':self.rad.tolist(),
            'bounds_ml':[.01,10.],'kernel':{'length_fraction':.3,'nugget':1e-8},'gas_ml_qp':1.,'bulge_ml_qp':.7,
            'optimized_direct_cache_qc':{'relative_to_peak':0.},
            'optimization':{'status':'feasible_weak_order','xi':[.5]*3,'optimality':{'passed':True},
                'max_scaled_violation':0.,'max_bound_violation':0.,'constraint_version':'positive_net_gbar_v1',
                'gbar_floor_relative':1e-8},'provenance':{'input_parser':columns.parser_provenance(),
                'script_sha256':inverse.sha256(inverse.__file__),'driver_sha256':inverse.sha256(inverse.mono.__file__),
                'utils_sha256':inverse.sha256(Path(inverse.DEFAULT_PROJECT_ROOT)/'utils.py')}}
        path=self.root/'corrected.json';path.write_text(json.dumps(data))
        xi,proof=gp.reference_profile(path,'Fixture',self.rad,.01,10.,args,{'fixture':'hash'})
        np.testing.assert_array_equal(xi,[.5]*3)
        data['provenance']['script_sha256']='legacy';path.write_text(json.dumps(data))
        with self.assertRaises(ValueError):gp.reference_profile(path,'Fixture',self.rad,.01,10.,args,{'fixture':'hash'})

    def test_plot_rejects_old_parent_and_draws_correct_decomposition(self):
        result={'galaxy':'Fixture','reference_ml':.5,'gas_ml_qp':1.,'bulge_ml_qp':.7,
            'retained_radii_kpc':[1.,2.],'optimization':{},'status':'infeasible','bounds_ml':[.01,20.],
            'input_sha256':{'rotmod':hashlib.sha256(self.rot.read_bytes()).hexdigest()},'provenance':{}}
        target=self.root/'plot.png'
        with self.assertRaisesRegex(ValueError,'legacy/unverified'):
            plotter.render_plot(result,target,self.root,self.root,inverse)
        self.assertFalse(target.exists())
        result['provenance']['input_parser']=columns.parser_provenance()
        with patch.object(inverse,'galaxy_data',return_value=self.load()),patch.object(plotter.plt,'close'):
            qc=plotter.render_plot(result,target,self.root,self.root,inverse)
            figure=plotter.plt.gcf();axis=figure.axes[0]
            disc=next(line for line in axis.lines if line.get_label().startswith('SPARC disc'))
            gas=next(line for line in axis.lines if line.get_label()=='Gas $g_{gas}$')
            np.testing.assert_allclose(disc.get_ydata(),inverse.vel2acc(.5*np.array([40.,50.])**2,np.array([1.,2.])))
            np.testing.assert_allclose(gas.get_ydata(),inverse.vel2acc(np.array([-16.,25.]),np.array([1.,2.])))
        plotter.plt.close(figure);self.assertTrue(target.exists())

if __name__=='__main__':unittest.main()

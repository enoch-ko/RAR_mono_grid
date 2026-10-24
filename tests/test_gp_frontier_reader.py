import hashlib,json,math,tempfile,unittest
from pathlib import Path
import numpy as np
from gp_frontier_reader import read_frontier_result,UnsupportedFrontierLayout
from plot_ml_profile_clouds import choose_record,load_gp,threshold_probability


class FrontierReaderTests(unittest.TestCase):
    def test_streamed_indices_hash_and_chunk_boundaries(self):
        records=[dict(valid=True,J_ref=float(i),marginal_logp=math.log(.01+i/1000),
                      source='SPARC_constant_0.5' if i==0 else 'GP_direction') for i in range(90)]
        data=dict(galaxy='synthetic',settings={},records=records,
                  frontier=[dict(record_index=i,xi=[.5,.6]) for i in [0,15,40,89]],provenance={})
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'result.json';path.write_text(json.dumps(data,indent=2)+'\n')
            for chunk in [17,127,1024*1024]:
                saved,digest=read_frontier_result(path,chunk)
                self.assertEqual(digest,hashlib.sha256(path.read_bytes()).hexdigest())
                self.assertEqual(saved['selection_reader']['evaluated_records'],90)
                self.assertEqual([r['original_candidate_index'] for r in saved['records']],[0,15,40,89])
                self.assertEqual([r['record_index'] for r in saved['frontier']],[0,1,2,3])
                self.assertEqual(choose_record(saved['records'],.05),(40.,2))

    def test_completion_count_mismatch_rejected(self):
        data=dict(galaxy='synthetic',records=[dict(valid=True,J_ref=0.,marginal_logp=-1.)],
                  frontier=[dict(record_index=0,xi=[.5])])
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)/'galaxies';folder.mkdir();path=folder/'synthetic_gp_search.json'
            path.write_text(json.dumps(data,indent=2)+'\n')
            (Path(tmp)/'synthetic_completion.json').write_text(json.dumps({'records':2}))
            with self.assertRaises(ValueError):read_frontier_result(path)

    def test_unrecognized_layout_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'result.json';path.write_text(json.dumps({'records':[], 'frontier':[]}))
            with self.assertRaises(UnsupportedFrontierLayout):read_frontier_result(path)

    def test_near_cut_fallback_recovers_skipped_candidate(self):
        cut=threshold_probability(2.,'one-sided')
        records=[dict(valid=True,J_ref=float(i),marginal_logp=math.log(p),source='GP_direction')
                 for i,p in enumerate([cut-1e-15,cut+1e-15])]
        data=dict(status='finite_search_completed',galaxy='synthetic',settings={},physical_domain={},
                  knots_kpc=[1.,2.],records=records,frontier=[dict(record_index=0,xi=[.5,.5])],
                  provenance=dict(source_sha256={},script_sha256='test',input_sha256={}))
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);ck=root/'synthetic_gp_checkpoints';ck.mkdir()
            profile_path=ck/'profiles.npy';np.save(profile_path,np.array([[.5,.5],[.6,.6]]))
            data['checkpoint_storage']=dict(mode='compact',shape=[2,2],sha256=hashlib.sha256(profile_path.read_bytes()).hexdigest())
            path=root/'synthetic_gp_search.json';path.write_text(json.dumps(data,indent=2)+'\n')
            panels,report=load_gp([root],[2.],'one-sided','frontier-stream')
            self.assertEqual(next(iter(panels.values()))[0]['J_ref'],1.)
            self.assertEqual(report['source_files'][0]['selection_reader']['mode'],'full_near_cut_fallback')


if __name__=='__main__':unittest.main()

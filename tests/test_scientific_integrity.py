"""Check information partitions, genomic coding and simulation accounting."""
from pathlib import Path
import json
import sys

import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from uavbreed.models import GenomicEncoder,fit_ridge
from uavbreed.programme import field_trial,split_masks,forecast,management_experiment,safe_corr

@pytest.fixture
def cohort():
    rng=np.random.default_rng(34)
    truth=pd.DataFrame({'genotype_id':[str(i) for i in range(12)],'family':[f'F{i//2}' for i in range(12)],'mother_id':[f'M{i//2}' for i in range(12)],'father_id':[f'P{i//2}' for i in range(12)]})
    for j in range(1,5):truth[f'gv{j}']=rng.normal(size=12)
    markers=pd.DataFrame(rng.binomial(2,.5,(12,14)),columns=[f'SNP{i}' for i in range(14)])
    markers.insert(0,'genotype_id',truth.genotype_id)
    cfg=json.loads((ROOT/'configs/pilot.json').read_text(encoding='utf-8'))
    crop='coffee' if 'coffee' in cfg['crops'] else 'cassava'
    trial,hidden=field_trial(truth,cfg,crop,98)
    return cfg,crop,truth,markers,trial,hidden

def test_mendelian_marker_support_and_disjoint_qtl_outputs():
    root=ROOT/'analysis/results/validated_pilot'
    files=list(root.glob('snapshots/**/genetics_audit.json'))
    if not files:pytest.skip('Run the programme first')
    for f in files:
        audit=json.loads(f.read_text(encoding='utf-8'))
        assert audit['ploidy']==2 and audit['qtl_marker_overlap']==0
        m=pd.read_csv(f.parent/'markers.csv').filter(regex='^SNP')
        assert set(np.unique(m))<=set([0,1,2])
        p=f.parent/'parent_markers.csv'
        if p.exists():assert set(np.unique(pd.read_csv(p).filter(regex='^SNP')))<=set([0,2])

def test_future_harvest_never_trains_on_test_year(cohort):
    _,crop,_,_,data,_=cohort
    for kind in ['future_harvest','new_parents','new_site','new_parents_site']:
        a,b=split_masks(data,crop,kind)
        assert not (a&b).any()
        assert data.loc[a,'year'].max()<data.loc[b,'year'].min()

def test_new_parent_holdout_excludes_both_sex_roles(cohort):
    _,crop,_,_,data,_=cohort
    a,b=split_masks(data,crop,'new_parents')
    parents_a=set(data.loc[a,'mother_id'])|set(data.loc[a,'father_id'])
    parents_b=set(data.loc[b,'mother_id'])|set(data.loc[b,'father_id'])
    assert not parents_a&parents_b

def test_test_targets_do_not_change_predictions_or_parent_index(cohort):
    _,crop,_,markers,data,_=cohort;a,b=split_masks(data,crop,'new_parents')
    pred,parent=forecast(data,markers,crop,a,b,'genomic_uav',return_parent=True)
    changed=data.copy();changed.loc[b,['yield_obs','quality_obs','disease_obs']]=12345
    other,other_parent=forecast(changed,markers,crop,a,b,'genomic_uav',return_parent=True)
    np.testing.assert_allclose(pred,other);np.testing.assert_allclose(parent,other_parent)

def test_parent_model_excludes_candidate_uav_direct_effects(cohort):
    _,crop,_,markers,data,_=cohort;a,b=split_masks(data,crop,'future_harvest')
    _,parent=forecast(data,markers,crop,a,b,'genomic_uav',return_parent=True)
    changed=data.copy();changed.loc[b,['height_early','cover_early','ndre_early','thermal_early','senescence_early']]=999
    _,other=forecast(changed,markers,crop,a,b,'genomic_uav',return_parent=True)
    np.testing.assert_allclose(parent,other)

def test_interpool_dominance_is_same_locus_product():
    x=np.array([[0,0],[0,2],[2,0],[2,2]],float)
    y=np.array([[0,2],[2,2],[0,0],[2,0]],float)
    enc=GenomicEncoder().fit(x,y);z=enc.transform(x,y)
    za=x/2-.5;zb=y/2-.5
    expected=-2*za*zb/np.sqrt(np.sum(4*np.full(2,.25**2)))
    np.testing.assert_allclose(z[:,-2:],expected)

def test_allele_frequencies_fit_training_only():
    enc=GenomicEncoder().fit(np.array([[0,1],[1,2],[2,0]]))
    before=enc.p.copy();enc.transform(np.array([[2,2]]))
    np.testing.assert_equal(enc.p,before)

def test_constant_prediction_has_undefined_rank_correlation():
    assert np.isnan(safe_corr(np.arange(5),np.ones(5)))

def test_field_seed_reproduces_observation_and_hidden_truth(cohort):
    cfg,crop,truth,_,data,hidden=cohort
    other,other_truth=field_trial(truth,cfg,crop,98)
    pd.testing.assert_frame_equal(data,other);pd.testing.assert_frame_equal(hidden,other_truth)

def test_cassava_harvest_uses_distinct_annual_plots(cohort):
    cfg,crop,_,_,data,_=cohort
    if crop!='cassava':pytest.skip('Cassava destructive harvest')
    assert data.groupby('plot_id').year.nunique().eq(1).all()
    counts=data.groupby(['genotype_id','site','management','rep']).plot_id.nunique()
    assert counts.eq(cfg['harvests']).all()

def test_oracle_bounds_management_policy_on_identical_fields(cohort):
    cfg,crop,*_=cohort
    records,_=management_experiment(cfg,crop,777)
    values={r['policy']:r['net_yield_equivalent'] for r in records}
    assert values['oracle']>=max(v for k,v in values.items() if k!='oracle')-1e-12

def test_saved_parent_holdouts_and_chronology():
    p=ROOT/'analysis/results/validated_pilot/split_audit.csv'
    if not p.exists():pytest.skip('Run the programme first')
    data=pd.read_csv(p)
    assert data.row_overlap.eq(0).all()
    assert (data.train_year_max<data.test_year_min).all()
    assert data.loc[data.validation.isin(['new_parents','new_parents_site']),'all_parent_overlap'].eq(0).all()

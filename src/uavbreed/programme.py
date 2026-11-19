"""Run paired breeding and randomized management experiments."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_squared_error,r2_score
from threadpoolctl import threadpool_limits

from .models import GenomicEncoder,fit_ridge

ENV=['water','n_supply','pressure','thermal_time','management']
UAV=['height_early','cover_early','ndre_early','thermal_early','senescence_early']
TARGETS=['yield_obs','quality_obs','disease_obs']

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def r_engine(root,config,crop,seed,state,folder,rlib,rscript,action):
    folder.mkdir(parents=True,exist_ok=True)
    cmd=[str(rscript),'--vanilla',str(root/'src/uavbreed/genetics.R'),action,str(root),str(config),crop,str(seed),str(state),str(folder),str(rlib)]
    result=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace')
    (folder/'engine.log').write_text(result.stdout+result.stderr,encoding='utf-8')
    if result.returncode:raise RuntimeError('Genetics engine failed: '+str(folder/'engine.log'))
    audit=json.loads((folder/'genetics_audit.json').read_text(encoding='utf-8'))
    if audit['qtl_marker_overlap']!=0:raise ValueError('Marker/QTL overlap')
    return pd.read_csv(folder/'genetic_truth.csv',dtype={'genotype_id':str,'mother_id':str,'father_id':str}),pd.read_csv(folder/'markers.csv',dtype={'genotype_id':str})

def field_trial(truth,cfg,crop,seed):
    rng=np.random.default_rng(seed);c=cfg['crops'][crop]
    n=len(truth)
    grid=pd.MultiIndex.from_product([range(n),range(cfg['sites']),range(cfg['harvests']),range(cfg['management_levels']),range(cfg['field_replicates'])],names=['g','site','year','management','rep']).to_frame(index=False)
    g=grid.g.to_numpy();site=grid.site.to_numpy();year=grid.year.to_numpy();m=grid.management.to_numpy()
    t=truth.iloc[g].reset_index(drop=True)
    for column in ['genotype_id','family','mother_id','father_id']:grid[column]=t[column].to_numpy()
    positions={}
    def position(row):
        planting_year=0 if crop=='coffee' else int(row.year)
        key=(int(row.site),int(row.management),int(row.rep),planting_year)
        if key not in positions:positions[key]=rng.permutation(n)
        return int(positions[key][int(row.g)])
    grid['plot_position']=[position(row) for row in grid.itertuples()]
    grid['plot_id']=[f'{crop}_S{r.site}_M{r.management}_B{r.rep}_G{r.g}_Y{0 if crop=="coffee" else r.year}' for r in grid.itertuples()]
    water=.25+.12*site+.10*year-.16*m
    nitrogen=.4+.2*m+.05*site
    pressure=.15+.08*site+.06*year
    time_sum=1+.1*year-.05*site
    grid['water']=water;grid['n_supply']=nitrogen;grid['pressure']=pressure;grid['thermal_time']=time_sum
    a=t.gv1.to_numpy();quality=t.gv2.to_numpy();response=t.gv3.to_numpy();resistance=t.gv4.to_numpy()
    genotype_environment=a+.30*response*water+.20*a*nitrogen
    disease=1/(1+np.exp(-(pressure*2-.65*resistance-1)))
    potential=c['mean_yield']+c['yield_scale']*(genotype_environment-.8*water+.3*nitrogen-.5*disease)
    if crop=='coffee':potential+=c['yield_scale']*.25*np.cos(year*np.pi+a)
    block=rng.normal(0,.12,(cfg['sites'],cfg['harvests'],cfg['management_levels'],cfg['field_replicates']))
    permanent=rng.normal(0,.12,(n,cfg['sites'],cfg['management_levels'],cfg['field_replicates'])) if crop=='coffee' else np.zeros((n,cfg['sites'],cfg['management_levels'],cfg['field_replicates']))
    noise=rng.normal(0,cfg['phenotyping']['residual_sd'],len(grid))
    grid['yield_obs']=potential+c['yield_scale']*(noise+block[site,year,m,grid.rep.to_numpy()]+permanent[g,site,m,grid.rep.to_numpy()])
    grid['quality_obs']=c['quality_mean']+c['quality_scale']*(quality+.10*water+rng.normal(0,.25,len(grid)))
    grid['disease_obs']=np.clip(disease+rng.normal(0,.05,len(grid)),0,1)
    sd=cfg['phenotyping']['sensor_noise_sd']
    grid['height_early']=1+.22*a+.12*quality-.15*water+rng.normal(0,sd,len(grid))
    grid['cover_early']=1/(1+np.exp(-(.35*a+.25*response-.4*water+rng.normal(0,sd,len(grid)))))
    grid['ndre_early']=.35+.08*a+.06*nitrogen-.06*disease+rng.normal(0,sd*.3,len(grid))
    grid['thermal_early']=3*water-.4*response+rng.normal(0,sd,len(grid))
    grid['senescence_early']=pressure+.2*quality+.5*disease+rng.normal(0,sd,len(grid))
    missing=rng.random((len(grid),len(UAV)))<cfg['phenotyping']['missing_probability']
    for j,col in enumerate(UAV):grid.loc[missing[:,j],col]=np.nan
    truth_rows=pd.DataFrame({'genotype_id':t.genotype_id,'potential_yield':potential,'potential_quality':c['quality_mean']+c['quality_scale']*quality,'potential_disease':disease})
    return grid,truth_rows

def genomic_arrays(trial,markers,parents=None):
    columns=[c for c in markers if c.startswith('SNP')]
    mat=markers.set_index('genotype_id')[columns]
    if parents is None:return mat.loc[trial.genotype_id].to_numpy(),None
    parent=parents.set_index('genotype_id')[columns]
    return parent.loc[trial.mother_id].to_numpy(),parent.loc[trial.father_id].to_numpy()

def split_masks(data,crop,kind):
    last_year=data.year.max();last_site=data.site.max()
    if crop=='maize':
        mothers=sorted(data.mother_id.unique());fathers=sorted(data.father_id.unique())
        hm=set(mothers[-max(2,len(mothers)//4):]);hf=set(fathers[-max(2,len(fathers)//4):])
        held_m=data.mother_id.isin(hm);held_f=data.father_id.isin(hf)
        new=held_m&held_f;unseen=held_m|held_f
    else:
        families=sorted(data.family.unique());held=set(families[-max(2,len(families)//4):])
        one=data.loc[data.family==families[-1],['mother_id','father_id']].iloc[0]
        held_parents=set(one.tolist())
        new=data.mother_id.isin(held_parents)&data.father_id.isin(held_parents)
        unseen=data.mother_id.isin(held_parents)|data.father_id.isin(held_parents)
    if kind=='future_harvest':train=data.year<last_year;test=data.year==last_year
    elif kind=='new_parents':train=(data.year<last_year)&~unseen;test=(data.year==last_year)&new
    elif kind=='new_site':train=(data.year<last_year)&(data.site<last_site);test=(data.year==last_year)&(data.site==last_site)
    elif kind=='new_parents_site':train=(data.year<last_year)&(data.site<last_site)&~unseen;test=(data.year==last_year)&(data.site==last_site)&new
    else:raise ValueError(kind)
    if (train&test).any() or not train.any() or not test.any():raise ValueError('Invalid split')
    return train.to_numpy(),test.to_numpy()

def forecast(data,markers,crop,train,test,kind,parents=None,return_parent=False):
    x1,x2=genomic_arrays(data,markers,parents)
    encoder=GenomicEncoder().fit(x1[train],None if x2 is None else x2[train])
    genomic=encoder.transform(x1,None if x2 is None else x2)
    additive=encoder.transform(x1,None if x2 is None else x2,additive_only=True)
    env=data[ENV].to_numpy();y=data[TARGETS].to_numpy()
    if kind=='field':
        means=data.loc[train].groupby('genotype_id')[TARGETS].mean()
        pred=data.loc[test,['genotype_id']].merge(means,on='genotype_id',how='left')[TARGETS].fillna(data.loc[train,TARGETS].mean()).to_numpy()
        parent_component=pred-data.loc[train,TARGETS].mean().to_numpy()
    else:
        z=np.c_[genomic,env]
        if kind=='genomic_uav':z=np.c_[z,data[UAV].to_numpy()]
        model=fit_ridge(z[train],y[train]);pred=model.predict(z[test])
        parent_component=None
        if return_parent:
            parent_model=fit_ridge(np.c_[additive,env][train],y[train])
            coefficients=parent_model.named_steps['ridge'].coef_[:,:additive.shape[1]]
            parent_component=additive[test]@coefficients.T
    return (pred,parent_component) if return_parent else pred

def safe_corr(a,b):
    return float(spearmanr(a,b).statistic) if np.std(a)>1e-12 and np.std(b)>1e-12 else float('nan')

def choose_parents(pred,trial,crop,cfg):
    scores=pred.groupby('genotype_id').parent_score.mean().sort_values(ascending=False)
    lookup=trial.drop_duplicates('genotype_id').set_index('genotype_id')
    if crop=='maize':
        candidates=lookup.loc[scores.index,['mother_id','father_id']].copy();candidates['score']=scores.to_numpy()
        a=candidates.groupby('mother_id').score.mean().nlargest(cfg['parent_count']//2).index.tolist()
        b=candidates.groupby('father_id').score.mean().nlargest(cfg['parent_count']//2).index.tolist()
        return a+b
    chosen=[];counts={}
    for ident in scores.index:
        fam=lookup.loc[ident,'family']
        if counts.get(fam,0)<2:
            chosen.append(ident);counts[fam]=counts.get(fam,0)+1
        if len(chosen)==cfg['parent_count']:break
    if len(chosen)!=cfg['parent_count']:raise ValueError('Parent family cap infeasible')
    return chosen

def management_experiment(cfg,crop,seed):
    rng=np.random.default_rng(seed);c=cfg['crops'][crop];mc=cfg['management']
    n=mc['training_fields']+mc['test_fields']
    stress=rng.uniform(0,1,n);response=rng.normal(0,1,n);sensor=stress+rng.normal(0,cfg['phenotyping']['sensor_noise_sd'],n)
    response_proxy=response+rng.normal(0,.4,n)
    x=np.c_[sensor,response_proxy,rng.normal(0,1,n)]
    effect=.9*stress+.15*response-.15
    effect[rng.random(n)<mc['nonresponse_rate']]=-.12
    y0=c['mean_yield']+c['yield_scale']*(.2*response-stress)
    y1=y0+c['yield_scale']*effect
    treatment=rng.integers(0,2,n)
    observed=np.where(treatment==1,y1,y0)+rng.normal(0,c['yield_scale']*.35,n)
    ntrain=mc['training_fields'];training=np.arange(n)<ntrain;testing=~training
    models=[]
    for action in [0,1]:
        mask=training&(treatment==action)
        model=HistGradientBoostingRegressor(max_iter=60,max_leaf_nodes=7,l2_regularization=2,random_state=seed).fit(x[mask],observed[mask]);models.append(model)
    cost=mc['treatment_cost_yield_units']*c['yield_scale']
    choice=models[1].predict(x[testing])-models[0].predict(x[testing])>cost
    policies={'none':np.zeros(testing.sum(),bool),'all':np.ones(testing.sum(),bool),'uav_rule':x[testing,0]>.5,'learned_policy':choice,'oracle':y1[testing]-y0[testing]>cost}
    records=[]
    for name,decision in policies.items():
        utility=np.where(decision,y1[testing],y0[testing])-cost*decision
        records.append({'crop':crop,'policy':name,'net_yield_equivalent':float(utility.mean()),'treated_fraction':float(decision.mean()),'fields':int(testing.sum()),'data_status':'SYNTHETIC'})
    return records,pd.DataFrame({'stress_truth':stress[testing],'sensor':sensor[testing],'response_proxy_observed':response_proxy[testing],'y0_truth':y0[testing],'y1_truth':y1[testing],'learned_action':choice})

def run(root,config_path,rscript,rlib,output=None):
    started=time.time();cfg=json.loads(config_path.read_text(encoding='utf-8'))
    out=output or root/'analysis/results/pilot';out.mkdir(parents=True,exist_ok=True)
    (out/'config_used.json').write_text(json.dumps(cfg,indent=2)+'\n',encoding='utf-8')
    metrics=[];selection=[];management=[];audits=[]
    for rep in range(cfg['replicates']):
        for ci,crop in enumerate(cfg['crops']):
            seed=cfg['seed']+rep*1000+ci*100
            base=out/f'snapshots/r{rep:03d}/{crop}';base.mkdir(parents=True,exist_ok=True)
            initial=base/'founders.rds'
            r_engine(root,config_path,crop,seed,initial,base/'initial',rlib,rscript,'init')
            for strategy in cfg['strategies']:
                state=base/(strategy+'.rds');shutil.copy2(initial,state)
                folder=base/'initial'
                for cycle in range(cfg['cycles']):
                    truth=pd.read_csv(folder/'genetic_truth.csv',dtype={'genotype_id':str,'mother_id':str,'father_id':str})
                    markers=pd.read_csv(folder/'markers.csv',dtype={'genotype_id':str})
                    parents=pd.read_csv(folder/'parent_markers.csv',dtype={'genotype_id':str}) if crop=='maize' else None
                    data,hidden=field_trial(truth,cfg,crop,seed+cycle*10+1)
                    for validation in ['future_harvest','new_parents','new_site','new_parents_site']:
                        train,test=split_masks(data,crop,validation)
                        returned=forecast(data,markers,crop,train,test,strategy,parents,return_parent=validation=='future_harvest')
                        pred,parent_component=returned if validation=='future_harvest' else (returned,None)
                        observed=data.loc[test,'yield_obs'].to_numpy()
                        metrics.append({'replicate':rep,'crop':crop,'strategy':strategy,'cycle':cycle,'validation':validation,'n_train':int(train.sum()),'n_test':int(test.sum()),'rmse':float(np.sqrt(mean_squared_error(observed,pred[:,0]))),'r2':float(r2_score(observed,pred[:,0])),'spearman':safe_corr(observed,pred[:,0])})
                        train_parents=set(data.loc[train,'mother_id'])|set(data.loc[train,'father_id'])
                        test_parents=set(data.loc[test,'mother_id'])|set(data.loc[test,'father_id'])
                        audits.append({'replicate':rep,'crop':crop,'strategy':strategy,'cycle':cycle,'validation':validation,'row_overlap':int((train&test).sum()),'genotype_overlap':len(set(data.loc[train,'genotype_id'])&set(data.loc[test,'genotype_id'])),'maternal_overlap':len(set(data.loc[train,'mother_id'])&set(data.loc[test,'mother_id'])),'paternal_overlap':len(set(data.loc[train,'father_id'])&set(data.loc[test,'father_id'])),'all_parent_overlap':len(train_parents&test_parents),'train_year_max':int(data.loc[train,'year'].max()),'test_year_min':int(data.loc[test,'year'].min())})
                        if validation=='future_harvest':
                            scales=data.loc[train,TARGETS].std().to_numpy();center=data.loc[train,TARGETS].mean().to_numpy();weights=np.array([.65,.20,-.15])
                            scores=((pred-center)/scales)@weights
                            selected=data.loc[test,['genotype_id','mother_id','father_id','family']].copy();selected['score']=scores
                            selected['parent_score']=(parent_component/scales)@weights
                            ids=selected.groupby('genotype_id').score.mean().nlargest(max(1,int(len(truth)*cfg['selection_fraction']))).index
                            h=hidden.loc[test].groupby('genotype_id').mean(numeric_only=True)
                            selection.append({'replicate':rep,'crop':crop,'strategy':strategy,'cycle':cycle,'selected_count':len(ids),'selected_potential_yield':float(h.loc[ids].potential_yield.mean()),'cohort_potential_yield':float(h.potential_yield.mean()),'mean_genetic_value':float(truth.gv1.mean()),'genetic_sd':float(truth.gv1.std()),'heterozygosity':float(((markers.filter(regex='^SNP')==1).to_numpy()).mean()),'calendar_year_assumption':cycle*cfg['crops'][crop]['cycle_years_assumed']})
                            parent_ids=choose_parents(selected,data,crop,cfg)
                            selected.to_csv(base/f'{strategy}_cycle{cycle}_scores.csv',index=False)
                            if crop=='coffee':
                                chosen=selected.groupby('genotype_id').score.mean().nlargest(6).index.tolist()
                                typed=truth.set_index('genotype_id')[['S1','S2']]
                                compat=[]
                                for recipient in chosen:
                                    female=set(typed.loc[recipient].tolist())
                                    for donor in chosen:
                                        accepted=float(np.mean([allele not in female for allele in typed.loc[donor]])) if donor!=recipient else 0.0
                                        compat.append({'recipient':recipient,'pollen_donor':donor,'synthetic_compatible_pollen_fraction':accepted,'planting_share_assumption':1/6,'status':'candidate set; compatibility proxy only, not measured fruit set'})
                                pd.DataFrame(compat).to_csv(base/f'{strategy}_cycle{cycle}_coffee_compatibility.csv',index=False)
                    if rep==0 and cycle==0 and strategy==cfg['strategies'][0]:
                        data.to_csv(base/'observations_SYNTHETIC.csv',index=False)
                        hidden.to_csv(base/'potential_outcomes_TRUTH_ONLY.csv',index=False)
                    if cycle+1<cfg['cycles']:
                        folder=base/f'{strategy}_cycle{cycle+1}';folder.mkdir(parents=True,exist_ok=True)
                        (folder/'parent_selection.json').write_text(json.dumps({'parent_ids':parent_ids}),encoding='utf-8')
                        r_engine(root,config_path,crop,seed+(cycle+1)*10,state,folder,rlib,rscript,'advance')
            records,raw=management_experiment(cfg,crop,seed+99)
            for record in records:record['replicate']=rep
            management.extend(records)
            if rep==0:raw.to_csv(base/'management_TRUTH_ONLY.csv',index=False)
            print(f'Completed replicate {rep+1}/{cfg["replicates"]}: {crop}',flush=True)
    for name,records in [('prediction_metrics',metrics),('selection_metrics',selection),('management_metrics',management),('split_audit',audits)]:
        pd.DataFrame(records).to_csv(out/(name+'.csv'),index=False)
    summary={'data_status':'SYNTHETIC_NOT_FIELD_OBSERVATIONS','replicates':cfg['replicates'],'cycles':cfg['cycles'],'crops':list(cfg['crops']),'seconds':round(time.time()-started,2),'config_sha256':digest(config_path),'engine_sha256':digest(root/'src/uavbreed/genetics.R'),'programme_sha256':digest(Path(__file__)),'interpretation':'Monte Carlo programme-design experiment; uncalibrated parameters, no commercial superiority or realised field genetic gain claimed.'}
    (out/'run_manifest.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    return out

from pathlib import Path
import os
import json,gzip,hashlib
import numpy as np
import pandas as pd
O=Path(os.environ.get("ANALYSIS_OUTPUT","outputs/analysis"))
PREP=Path(os.environ.get("PREP_OUTPUT","outputs/preparation"))
def save(x,p,index=False):x.to_csv(O/p,sep="\t",index=index)
spec=pd.read_csv(O/"cohort_specifications.tsv",sep="\t")
primary="primary7_min10_no_wu16"; original="qc_original8_min10_with_wu16"
res={r.analysis_id:pd.read_csv(O/r.analysis_id/"all_genes_results.tsv",sep="\t").set_index("ensembl_id",drop=False) for r in spec.itertuples(index=False)}
p=res[primary].copy();orig=res[original]
# Attach primary per-patient descriptive expression differences and cellular detection summaries.
pair=pd.read_csv(O/primary/"patient_logCPM_differences.tsv",sep="\t",index_col=0)
for patient in pair:p["paired_logCPM_delta_"+patient]=pair[patient]
p["n_patients_direction_agree"]=(pair.apply(np.sign).eq(np.sign(p.loc[pair.index,"logFC"]),axis=0)).sum(axis=1)
p["median_patient_logCPM_delta"]=pair.median(axis=1)
cohort_samples=spec.set_index("analysis_id").loc[primary,"samples"].split(";")
q=pd.read_csv(PREP/"retained_epithelial_cells.tsv.gz",sep="\t")
f=pd.read_csv(PREP/"feature_mapping.tsv",sep="\t")
ng=len(f);det={};den={}
for sample in cohort_samples:
 cells=q[q.sample_id.eq(sample)]
 lookup=np.full(737280,-1,dtype=np.int64)
 status=cells.epithelial_status.map({"Normal":0,"Malignant":1}).to_numpy()
 lookup[cells.raw_matrix_column_1based.to_numpy()-1]=status
 accum=np.zeros((ng,2),dtype=np.int64)
 file=next(Path("data").glob("*_"+sample+"_matrix.mtx.gz"))
 with gzip.open(file,"rt") as stream:
  next(stream);line=next(stream)
  while line.startswith("%") :line=next(stream)
  for c in pd.read_csv(stream,sep=r"\s+",header=None,names=["row","col","value"],dtype=np.int64,chunksize=1000000):
   group=lookup[c.col.to_numpy()-1];keep=(group>=0)&(c.value.to_numpy()>0)
   ix=(c.row.to_numpy()[keep]-1)*2+group[keep]
   accum+=np.bincount(ix,minlength=ng*2).reshape(ng,2)
 for j,name in enumerate(["Normal","Malignant"]):
  key=sample+"__"+name;det[key]=accum[:,j];den[key]=int((status==j).sum())
detect=pd.DataFrame(det,index=f.ensembl_id);detect.index.name="ensembl_id"
save(detect,"primary7_cellular_detection_counts.tsv.gz",True)
save(pd.DataFrame({"pseudobulk_id":list(den),"n_cells":list(den.values())}),"primary7_cellular_detection_denominators.tsv")
for name in ["Normal","Malignant"]:
 cols=[c for c in detect if c.endswith("__"+name)]
 p["cell_detection_fraction_"+name]=detect[cols].sum(axis=1)/sum(den[c] for c in cols)
 p["patients_detected_"+name]=(detect[cols]>0).sum(axis=1)
# LOO fits are effect-only and have no new FDR family.
influence_summary=[]
for parent,prefix,label in [(primary,"loo_primary7_","primary7"),(original,"loo_original8_","original8")]:
 base=res[parent];base=base[base.logFC.notna()]
 ids=[x for x in res if x.startswith(prefix)]
 mat=pd.DataFrame({x:res[x].logFC for x in ids}).reindex(base.index)
 stab=pd.DataFrame(index=base.index)
 stab["ensembl_id"]=base.index;stab["gene_symbol"]=base.gene_symbol;stab["parent_logFC"]=base.logFC;stab["parent_FDR"]=base.FDR
 stab["n_loo_estimable"]=mat.notna().sum(axis=1)
 stab["n_loo_same_direction"]=np.sign(mat).eq(np.sign(base.logFC),axis=0).sum(axis=1)
 stab["all_loo_same_direction"]=(stab.n_loo_same_direction==len(ids))
 stab["loo_logFC_min"]=mat.min(axis=1);stab["loo_logFC_max"]=mat.max(axis=1)
 delta=mat.subtract(base.logFC,axis=0).abs()
 stab["max_abs_loo_effect_change"]=delta.max(axis=1);stab["most_influential_omission"]=delta.idxmax(axis=1)
 stab=pd.concat([stab,mat],axis=1)
 save(stab,"influence_"+label+"_gene_summary.tsv")
 if parent==primary:
  for c in ["n_loo_estimable","n_loo_same_direction","all_loo_same_direction","loo_logFC_min","loo_logFC_max","max_abs_loo_effect_change","most_influential_omission"]:p[c]=stab[c]
 for name in ids:
  valid=base.logFC.notna()&res[name].logFC.reindex(base.index).notna()
  sig=base.FDR.lt(.05)
  influence_summary.append(dict(parent=parent,omitted_fit=name,n_genes_compared=int(valid.sum()),effect_Pearson=base.logFC[valid].corr(res[name].logFC.reindex(base.index)[valid]),median_abs_effect_change=float(delta[name].median()),max_abs_effect_change=float(delta[name].max()),parent_significant_genes=int(sig.sum()),sign_flips_in_parent_significant=int((sig & np.sign(base.logFC).ne(np.sign(mat[name]))).sum())))
save(pd.DataFrame(influence_summary),"influence_cohort_summary.tsv")
comparison=[]
for name in spec[spec.role.isin(["qc_sensitivity","prespecified_sensitivity"])].analysis_id:
 s=res[name];valid=p.logFC.notna()&s.logFC.notna();sig=p.FDR.lt(.05);co=s.FDR.lt(.05)
 p[name+"__logFC"]=s.logFC;p[name+"__FDR"]=s.FDR;p[name+"__test_status"]=s.test_status
 comparison.append(dict(analysis_id=name,n_common_tested=int(valid.sum()),effect_Pearson=p.loc[valid,"logFC"].corr(s.loc[valid,"logFC"]),primary_sig_tested=int((sig&valid).sum()),primary_sig_also_FDR05=int((sig&co).sum()),primary_sig_same_direction=int((sig&valid&np.sign(p.logFC).eq(np.sign(s.logFC))).sum()),primary_sig_opposite_direction_FDR05=int((sig&co&np.sign(p.logFC).ne(np.sign(s.logFC))).sum())))
save(pd.DataFrame(comparison),"sensitivity_vs_primary.tsv")
save(p,"primary_all_genes_annotated.tsv")
sig=p[p.FDR.lt(.05)].sort_values(["FDR","logFC"],ascending=[True,False])
save(sig,"primary_significant_annotated.tsv")
save(sig[sig.abs_log2FC_ge_1],"primary_prioritized_FDR05_absFC1.tsv")
# Presentation examples selected deterministically from primary FDR only, not sensitivities.
examples=pd.concat([sig[sig.logFC>0].head(8),sig[sig.logFC<0].head(8)])
save(examples,"primary_example_genes.tsv")
for field in ["FDR","logFC"]:assert p[field].equals(res[primary][field])
assert len(p)==33538 and p.logFC.notna().sum()==6720
summary={"primary_tested":int(p.logFC.notna().sum()),"primary_excluded_abundance":int(p.logFC.isna().sum()),"primary_FDR05":len(sig),"primary_up":int((sig.logFC>0).sum()),"primary_down":int((sig.logFC<0).sum()),"primary_FDR05_absFC1":int(sig.abs_log2FC_ge_1.sum()),"sig_all7_loo_same_direction":int(sig.all_loo_same_direction.sum()),"sig_all7_patients_same_direction":int(sig.n_patients_direction_agree.eq(7).sum()),"sig_6plus_patients_same_direction":int(sig.n_patients_direction_agree.ge(6).sum()),"top_up_examples":examples[examples.logFC>0][["gene_symbol","logFC","FDR"]].to_dict("records"),"top_down_examples":examples[examples.logFC<0][["gene_symbol","logFC","FDR"]].to_dict("records")}
(O/"results_summary.json").write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2));print("\nSENSITIVITIES\n",pd.DataFrame(comparison).to_string(index=False));print("\nINFLUENCE\n",pd.DataFrame(influence_summary).to_string(index=False))

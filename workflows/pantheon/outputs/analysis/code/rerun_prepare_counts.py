import os
import gzip, hashlib, json, sys, platform, time
from pathlib import Path
import importlib.metadata as im
import numpy as np
import pandas as pd
from scipy import sparse

O=Path(os.environ["PREP_OUTPUT"])
O.mkdir(exist_ok=True,parents=True)
def save(df,name,index=False):
    df.to_csv(O/name,sep="\t",index=index)
def say(*x): print(*x,flush=True)
versions={"python":sys.version,"platform":platform.platform(),"packages":{d.metadata["Name"]:d.version for d in im.distributions()}}
(O/"python_versions.json").write_text(json.dumps(versions,indent=2))
# Input and previous-stage provenance; do not modify earlier artifacts.
manifest=[]
for base in [Path("data"),Path("outputs/inspection"),Path("outputs/planning_revision")]:
 for f in sorted(base.rglob("*")):
  if f.is_file():
   h=hashlib.sha256()
   with f.open("rb") as z:
    for chunk in iter(lambda:z.read(1<<20),b""):h.update(chunk)
   manifest.append(dict(path=str(f),bytes=f.stat().st_size,sha256=h.hexdigest()))
save(pd.DataFrame(manifest),"input_and_prior_artifact_checksums.tsv")
a=pd.read_csv("data/GSE242230_annotations.txt.gz",sep="\t")
assert not a.cell_id.duplicated().any()
meta=pd.read_csv(Path(os.environ["INSPECTION_OUTPUT"])/"joined_sample_metadata.tsv",sep="\t")
assert meta.sample_id.is_unique and meta.table_patient_id.is_unique
assert set(meta.sample_id)==set(a.sample_id)
meta=meta.rename(columns={"table_patient_id":"patient_id"})
meta["mapping_basis"]="Numeric WU patient identifier, supported by matching SOFT/Table S1 clinical fields"
for x,y in [("ajcc stage","table_stage"),("tumor location","table_location"),("pathology","table_pathology")]:
 assert meta[x].str.strip().equals(meta[y].str.strip())
save(meta,"sample_patient_metadata.tsv")
labels={"Malignant - Basal":"Malignant","Malignant - Classical":"Malignant","Normal Epithelial":"Normal"}
Q=[]; checks=[]; pb={}; subpb={}; pbmeta=[]; submeta=[]
features=None
for matrix in sorted(Path("data").glob("*_matrix.mtx.gz")):
 start=time.monotonic()
 prefix=str(matrix).removesuffix("_matrix.mtx.gz")
 sample=matrix.name.split("_")[1]
 say("START",sample)
 f=pd.read_csv(prefix+"_features.tsv.gz",sep="\t",header=None,names=["ensembl_id","gene_symbol","feature_type"])
 if features is None:
  features=f.copy()
  assert features.ensembl_id.is_unique
  features["matrix_row_1based"]=np.arange(1,len(features)+1)
  features["is_mitochondrial"]=features.gene_symbol.str.startswith("MT-")
  features["is_ribosomal_symbol"]=features.gene_symbol.str.match(r"^RP[SL][0-9]")
  save(features,"feature_mapping.tsv")
 else: assert f.equals(features[["ensembl_id","gene_symbol","feature_type"]])
 b=pd.read_csv(prefix+"_barcodes.tsv.gz",header=None,names=["barcode"])
 assert b.barcode.is_unique
 ann=a[a.sample_id.eq(sample)].copy().reset_index(drop=True)
 ann["barcode"]=ann.cell_id.str.removeprefix(sample+"_")
 assert (sample+"_"+ann.barcode).equals(ann.cell_id)
 cols=pd.Index(b.barcode).get_indexer(ann.barcode)
 assert (cols>=0).all() and len(np.unique(cols))==len(cols)
 ann["raw_matrix_column_1based"]=cols+1
 ann["patient_id"]=meta.set_index("sample_id").loc[sample,"patient_id"]
 lookup=np.full(len(b),-1,dtype=np.int64);lookup[cols]=np.arange(len(ann))
 with gzip.open(matrix,"rt") as stream:
  header=stream.readline().strip(); assert header=="%%MatrixMarket matrix coordinate integer general"
  line=stream.readline()
  while line.startswith("%") :line=stream.readline()
  ng,nc,nnz=map(int,line.split());assert ng==len(f) and nc==len(b)
  rr=[];cc=[];vv=[];entries=0;rawsum=0;zeros=0;maxcount=0
  for chunk in pd.read_csv(stream,sep=r"\s+",header=None,names=["row","col","value"],dtype=np.int64,chunksize=1000000):
   # A noninteger token or missing value fails parsing rather than being rounded.
   z=chunk.to_numpy();r=z[:,0]-1;c=z[:,1]-1;v=z[:,2]
   assert ((r>=0)&(r<ng)&(c>=0)&(c<nc)&(v>=0)).all()
   entries+=len(z);rawsum+=int(v.sum());zeros+=int((v==0).sum());maxcount=max(maxcount,int(v.max()))
   dest=lookup[c];keep=dest>=0
   rr.append(r[keep]);cc.append(dest[keep]);vv.append(v[keep])
  assert entries==nnz
  r=np.concatenate(rr);c=np.concatenate(cc);v=np.concatenate(vv)
  selected_nonzero=int(np.count_nonzero(v))
  x=sparse.coo_matrix((v,(r,c)),shape=(ng,len(ann)),dtype=np.int64).tocsc()
  x.sum_duplicates();x.eliminate_zeros()
  duplicate_nonzero_entries=selected_nonzero-x.nnz
  assert duplicate_nonzero_entries==0, "Duplicated annotated gene-cell coordinates require review"
  totals=np.asarray(x.sum(axis=0)).ravel();detected=x.getnnz(axis=0)
  mt=np.asarray(x[features.is_mitochondrial.to_numpy()].sum(axis=0)).ravel()
  # Gene dominance is descriptive QC, not DE-based feature selection.
  top1=np.zeros(len(ann),dtype=np.int64);top10=top1.copy();toprow=np.full(len(ann),-1,dtype=np.int64)
  for j in range(len(ann)):
   lo,hi=x.indptr[j:j+2];vals=x.data[lo:hi]
   if len(vals):
    k=int(np.argmax(vals));top1[j]=vals[k];toprow[j]=x.indices[lo+k]
    top10[j]=np.partition(vals,max(0,len(vals)-10))[-10:].sum()
  ann["top_gene_ensembl_id"]=[features.ensembl_id.iloc[k] if k>=0 else "" for k in toprow]
  ann["top_gene_symbol"]=[features.gene_symbol.iloc[k] if k>=0 else "" for k in toprow]
  ann["top_gene_fraction"]=np.divide(top1,totals,out=np.full(len(totals),np.nan),where=totals>0)
  ann["top10_genes_fraction"]=np.divide(top10,totals,out=np.full(len(totals),np.nan),where=totals>0)
  ann["total_umi"]=totals;ann["detected_genes"]=detected;ann["mitochondrial_umi"]=mt
  ann["mitochondrial_fraction"]=np.divide(mt,totals,out=np.full(len(totals),np.nan),where=totals>0)
  ann["epithelial_status"]=ann.cell_type_specific.map(labels)
  ann["selected_by_annotation"]=(~ann.filtered)&ann.epithelial_status.notna()
  ann["retained_for_pseudobulk"]=ann.selected_by_annotation & (totals>0)
  ann["exclusion_reason"]=np.select([ann.filtered,ann.epithelial_status.isna(),totals==0],["supplied_filtered_True","not_target_epithelial_label","zero_count_cell"],default="retained")
  # Descriptive review flags ONLY, not exclusion thresholds or a doublet classifier.
  ann["flag_lt200_genes"]=ann.detected_genes<200
  ann["flag_mt_gt20pct"]=ann.mitochondrial_fraction>0.20
  ann["flag_mt_gt50pct"]=ann.mitochondrial_fraction>0.50
  ann["flag_gt50000_umi"]=ann.total_umi>50000
  ann["flag_gt7000_genes"]=ann.detected_genes>7000
  Q.append(ann)
  for colname,outdict,outmeta in [("epithelial_status",pb,pbmeta),("cell_type_specific",subpb,submeta)]:
   for group in sorted(ann.loc[ann.retained_for_pseudobulk,colname].unique()):
    mask=(ann.retained_for_pseudobulk & ann[colname].eq(group)).to_numpy()
    ident=sample+"__"+group.replace(" - ","_").replace(" ","_")
    counts=np.asarray(x[:,mask].sum(axis=1)).ravel()
    assert counts.dtype==np.int64 and counts.sum()==totals[mask].sum()
    outdict[ident]=counts
    outmeta.append(dict(pseudobulk_id=ident,sample_id=sample,patient_id=ann.patient_id.iloc[0],group=group,n_cells=int(mask.sum()),raw_library_umi=int(counts.sum()),genes_with_counts=int(np.count_nonzero(counts)),mitochondrial_fraction=float(counts[features.is_mitochondrial].sum()/counts.sum())))
 checks.append(dict(sample_id=sample,n_features=ng,n_raw_barcodes=nc,n_entries=entries,total_raw_umi=rawsum,n_zero_entries=zeros,max_entry=maxcount,n_annotated=len(ann),n_matched=len(cols),annotated_nonzero_entries=x.nnz,duplicate_annotated_nonzero_entries=duplicate_nonzero_entries,annotated_umi=int(totals.sum()),selected_cells=int(ann.selected_by_annotation.sum()),retained_cells=int(ann.retained_for_pseudobulk.sum()),zero_count_selected=int((ann.selected_by_annotation & ann.total_umi.eq(0)).sum()),elapsed_seconds=time.monotonic()-start))
 say("COMPLETE",json.dumps(checks[-1]))
q=pd.concat(Q,ignore_index=True)
assert len(q)==len(a) and q.cell_id.is_unique
save(q,"all_annotated_cell_qc.tsv.gz")
ret=q[q.retained_for_pseudobulk].copy()
ret["pseudobulk_id"]=ret.sample_id+"__"+ret.epithelial_status
ret["subtype_pseudobulk_id"]=ret.sample_id+"__"+ret.cell_type_specific.str.replace(" - ","_",regex=False).str.replace(" ","_",regex=False)
save(ret,"retained_epithelial_cells.tsv.gz")
save(q[q.selected_by_annotation & ~q.retained_for_pseudobulk],"excluded_selected_cells.tsv")
save(pd.DataFrame(checks),"raw_matrix_validation.tsv")
pbm=pd.DataFrame(pbmeta).merge(meta,on=["sample_id","patient_id"],validate="many_to_one")
subm=pd.DataFrame(submeta).merge(meta,on=["sample_id","patient_id"],validate="many_to_one")
counts=pd.DataFrame(pb,index=features.ensembl_id);counts.index.name="ensembl_id"
scounts=pd.DataFrame(subpb,index=features.ensembl_id);scounts.index.name="ensembl_id"
assert list(counts.columns)==pbm.pseudobulk_id.tolist()
assert ret.total_umi.sum()==counts.to_numpy().sum()
assert ret.total_umi.sum()==scounts.to_numpy().sum()
save(counts,"raw_pseudobulk_counts_all.tsv.gz",True);save(scounts,"raw_subtype_pseudobulk_counts_all.tsv.gz",True)
save(pbm,"pseudobulk_metadata_all.tsv");save(subm,"subtype_pseudobulk_metadata_all.tsv")
# Summary for ALL annotated cells and each specific label (including filtered cells).
metrics=["total_umi","detected_genes","mitochondrial_fraction"]
summary=q.groupby(["sample_id","cell_type_specific","filtered"],dropna=False)[metrics].agg(["count","min",lambda x:x.quantile(.05),"median",lambda x:x.quantile(.95),"max"])
summary.columns=[m+"_"+stat.replace("<lambda_0>","q05").replace("<lambda_1>","q95") for m,stat in summary.columns]
save(summary,"cell_qc_by_sample_specific_label.tsv",True)
flags=[c for c in q if c.startswith("flag_")]
save(q.groupby(["sample_id","epithelial_status","filtered"],dropna=False)[flags].sum(),"qc_review_flag_counts.tsv",True)
# Cohort manifests use observed pseudobulks only: absent normal groups have no count column.
size=ret.groupby(["sample_id","epithelial_status"]).size().unstack(fill_value=0).reindex(meta.sample_id,fill_value=0)
subsize=ret.groupby(["sample_id","cell_type_specific"]).size().unstack(fill_value=0).reindex(meta.sample_id,fill_value=0)
elig=meta[["sample_id","patient_id","pathology"]].copy().set_index("sample_id")
elig["n_Malignant"]=size.Malignant;elig["n_Normal"]=size.Normal
elig["n_Classical"]=subsize["Malignant - Classical"];elig["n_Basal"]=subsize["Malignant - Basal"]
for n in [5,10,20]:
 elig["eligible_min"+str(n)]=(elig.n_Malignant>=n)&(elig.n_Normal>=n)
 elig["reason_min"+str(n)]=["eligible" if m>=n and norm>=n else ";".join((["Malignant_below_"+str(n)] if m<n else [])+(["Normal_below_"+str(n)] if norm<n else [])) for m,norm in zip(elig.n_Malignant,elig.n_Normal)]
elig["eligible_classical_min10"]=(elig.n_Classical>=10)&(elig.n_Normal>=10)
save(elig,"sample_eligibility.tsv",True)
cohorts={"primary_min10":elig.index[elig.eligible_min10].tolist(),"sensitivity_min5":elig.index[elig.eligible_min5].tolist(),"sensitivity_min20":elig.index[elig.eligible_min20].tolist(),"classical_min10":elig.index[elig.eligible_classical_min10].tolist()}
cohorts["sensitivity_min5_without_wu39"]=[s for s in cohorts["sensitivity_min5"] if s!="wu39"]
cohorts["sensitivity_min5_without_wu07"]=[s for s in cohorts["sensitivity_min5"] if s!="wu07"]
for sample in cohorts["primary_min10"]:cohorts["leave_out_"+sample]=[s for s in cohorts["primary_min10"] if s!=sample]
(O/"cohort_definitions.json").write_text(json.dumps(cohorts,indent=2))
manifest=[]
for name,samples in cohorts.items():
 classical=name=="classical_min10"
 cm=(subm if classical else pbm).copy()
 if classical:
  cm=cm[cm.group.isin(["Normal Epithelial","Malignant - Classical"])].copy()
  cm["status"]=cm.group.map({"Normal Epithelial":"Normal","Malignant - Classical":"Malignant"})
 else: cm["status"]=cm.group
 cm=cm[cm.sample_id.isin(samples)].copy()
 cm["status_code"]=cm.status.map({"Normal":0,"Malignant":1})
 cm=cm.sort_values(["sample_id","status_code"]).reset_index(drop=True)
 assert len(cm)==2*len(samples) and (cm.groupby("sample_id").status.nunique()==2).all()
 for _,r in cm.iterrows():manifest.append(dict(cohort=name,**r.to_dict()))
 if name in ["primary_min10","classical_min10"]:
  save(cm,name+"_metadata.tsv")
  cts=scounts if classical else counts
  save(cts.loc[:,cm.pseudobulk_id],name+"_raw_counts.tsv.gz",True)
save(pd.DataFrame(manifest),"cohort_pseudobulk_manifest.tsv")
primary=ret[ret.sample_id.isin(cohorts["primary_min10"])]
result={"scope":"Raw count validation, descriptive QC and pseudobulk preparation only. No gene filtering, normalization, dispersion estimation, fitting, DE tests or gene ranking.","n_annotated":len(q),"n_selected":int(q.selected_by_annotation.sum()),"n_retained_epithelial":len(ret),"n_zero_count_selected":int((q.selected_by_annotation&q.total_umi.eq(0)).sum()),"n_features":len(features),"n_observed_status_pseudobulks":len(pbm),"n_observed_subtype_pseudobulks":len(subm),"primary_samples":cohorts["primary_min10"],"primary_status_counts":primary.epithelial_status.value_counts().to_dict(),"primary_subtype_counts":primary.cell_type_specific.value_counts().to_dict(),"review_flags_not_exclusions":{"lt200_genes":"detected_genes < 200","mt_gt20pct":"mitochondrial_fraction > 0.20","mt_gt50pct":"mitochondrial_fraction > 0.50","gt50000_umi":"total_umi > 50000","gt7000_genes":"detected_genes > 7000"},"retained_flag_counts":ret[flags].sum().to_dict(),"primary_flag_counts":primary[flags].sum().to_dict(),"contrast":"statusMalignant in ~patient+status, Normal reference; positive log2FC is Malignant minus Normal"}
(O/"preparation_summary.json").write_text(json.dumps(result,indent=2))
say(json.dumps(result,indent=2));say("DONE; no DE inference executed.")

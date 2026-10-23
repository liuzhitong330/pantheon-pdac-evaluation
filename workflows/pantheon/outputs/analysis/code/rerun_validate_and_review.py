from pathlib import Path
import os
import gzip,json,hashlib
import numpy as np
import pandas as pd
O=Path(os.environ["PREP_OUTPUT"])
def save(x,name,index=False):x.to_csv(O/name,sep="\t",index=index)
q=pd.read_csv(O/"retained_epithelial_cells.tsv.gz",sep="\t")
f=pd.read_csv(O/"feature_mapping.tsv",sep="\t")
c=pd.read_csv(O/"raw_pseudobulk_counts_all.tsv.gz",sep="\t",index_col=0)
s=pd.read_csv(O/"raw_subtype_pseudobulk_counts_all.tsv.gz",sep="\t",index_col=0)
m=pd.read_csv(O/"pseudobulk_metadata_all.tsv",sep="\t")
assert c.index.equals(pd.Index(f.ensembl_id,name="ensembl_id"))
assert list(c)==m.pseudobulk_id.tolist() and np.issubdtype(c.to_numpy().dtype,np.integer)
assert (c.to_numpy()>=0).all() and (s.to_numpy()>=0).all()
assert np.array_equal(c.sum().values,m.raw_library_umi.values)
for _,row in m.iterrows():
 if row.group=="Malignant": names=[x for x in s if x.startswith(row.sample_id+"__Malignant_")]
 else:names=[row.sample_id+"__Normal_Epithelial"]
 assert np.array_equal(c[row.pseudobulk_id],s[names].sum(axis=1))
 # Exact retained identities partition each pseudobulk, and cell UMI totals match its counts.
 cells=q[q.pseudobulk_id.eq(row.pseudobulk_id)]
 assert len(cells)==row.n_cells and cells.total_umi.sum()==row.raw_library_umi
cohort=pd.read_csv(O/"cohort_pseudobulk_manifest.tsv",sep="\t")
designinfo=[]
for name,z in cohort.groupby("cohort",sort=False):
 z=z.copy();patients=z.patient_id.drop_duplicates().tolist()
 design=pd.DataFrame({"Intercept":np.ones(len(z),dtype=int)},index=z.pseudobulk_id)
 for pat in patients[1:]: design["patient"+pat]=(z.patient_id==pat).to_numpy().astype(int)
 design["statusMalignant"]=(z.status=="Malignant").to_numpy().astype(int)
 assert np.linalg.matrix_rank(design)==design.shape[1]
 assert (z.groupby("patient_id").status.nunique()==2).all()
 contrast=pd.Series(0,index=design.columns,name="Malignant_minus_Normal");contrast.loc["statusMalignant"]=1
 # For every patient, design row Malignant minus row Normal equals the contrast vector.
 for pat in patients:
  rows=z[z.patient_id==pat].set_index("status")
  delta=design.loc[rows.loc["Malignant","pseudobulk_id"]]-design.loc[rows.loc["Normal","pseudobulk_id"]]
  assert np.array_equal(delta.values,contrast.values)
 save(design,"design_"+name+".tsv",True)
 save(contrast.to_frame(),"contrast_"+name+".tsv",True)
 designinfo.append(dict(cohort=name,n_patients=len(patients),n_pseudobulks=len(z),rank=design.shape[1],nominal_residual_df=len(z)-design.shape[1],patient_reference=patients[0],status_reference="Normal",coefficient="statusMalignant",positive_effect="Malignant higher than Normal",patient_factor_levels=";".join(patients)))
save(pd.DataFrame(designinfo),"design_validation.tsv")
# Hemoglobin dominance is a count-composition QC finding, not a DE-selected gene set.
hbs=["HBA1","HBA2","HBB"]
hbids=f.loc[f.gene_symbol.isin(hbs),"ensembl_id"].tolist()
review=m[["pseudobulk_id","sample_id","patient_id","group","n_cells","raw_library_umi","mitochondrial_fraction"]].copy()
for symbol in hbs:
 ids=f.loc[f.gene_symbol.eq(symbol),"ensembl_id"]
 review[symbol+"_umi"]=c.loc[ids].sum().values
 review[symbol+"_fraction"]=review[symbol+"_umi"]/review.raw_library_umi
review["HBA1_HBA2_HBB_fraction"]=c.loc[hbids].sum().values/review.raw_library_umi
save(review,"hemoglobin_composition_all_pseudobulks.tsv")
# Audit wu16 directly against raw matrix to quantify cell-level HBB/HBA contributions.
w=q[q.sample_id.eq("wu16")].copy();positions=w.raw_matrix_column_1based.to_numpy()
lookup={int(col):j for j,col in enumerate(positions)}
hbrows={int(r.matrix_row_1based):r.gene_symbol for _,r in f[f.gene_symbol.isin(hbs)].iterrows()}
values={symbol:np.zeros(len(w),dtype=np.int64) for symbol in hbs}
p=next(Path("data").glob("*_wu16_matrix.mtx.gz"))
with gzip.open(p,"rt") as stream:
 next(stream);line=next(stream)
 while line.startswith("%") :line=next(stream)
 for chunk in pd.read_csv(stream,sep=r"\s+",header=None,names=["row","col","value"],dtype=np.int64,chunksize=1000000):
  relevant=chunk[chunk.row.isin(hbrows) & chunk.col.isin(lookup)]
  for r in relevant.itertuples(index=False):values[hbrows[r.row]][lookup[r.col]]+=r.value
for sym,v in values.items():w[sym+"_umi"]=v;w[sym+"_fraction"]=v/w.total_umi
w["HBA1_HBA2_HBB_fraction"]=sum(values.values())/w.total_umi
w["review_flag_HBB_gt50pct"]=w.HBB_fraction>0.5
save(w,"wu16_cell_composition_review.tsv.gz")
flags=q.groupby(["sample_id","epithelial_status"]).agg(n_cells=("cell_id","size"),median_umi=("total_umi","median"),median_genes=("detected_genes","median"),median_top_gene_fraction=("top_gene_fraction","median"),median_top10_fraction=("top10_genes_fraction","median"),median_mt_fraction=("mitochondrial_fraction","median"),n_under500_genes=("detected_genes",lambda x:int((x<500).sum())),n_topgene_gt50pct=("top_gene_fraction",lambda x:int((x>.5).sum())))
save(flags,"selected_group_composition_qc.tsv",True)
wsummary=w.groupby("epithelial_status").agg(n_cells=("cell_id","size"),median_umi=("total_umi","median"),median_genes=("detected_genes","median"),median_HBB_fraction=("HBB_fraction","median"),n_HBB_gt50pct=("review_flag_HBB_gt50pct","sum"),median_HBA_HBB_fraction=("HBA1_HBA2_HBB_fraction","median"))
save(wsummary,"wu16_composition_summary.tsv",True)
print("Independent saved-count, subtype-sum and metadata assertions passed.")
print(pd.DataFrame(designinfo).to_string(index=False))
print("\nWU16 CELL QC\n"+wsummary.to_string())
print("\nWU16 PSEUDOBULK COMPOSITION\n"+review[review.sample_id.eq("wu16")].to_string(index=False))
print("\nGROUP QC\n"+flags.to_string())
print("No DE tests, normalization or gene selection performed. No cells excluded for these review flags.")

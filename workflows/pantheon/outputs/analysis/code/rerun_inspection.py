from pathlib import Path
import os
import gzip, hashlib, json, platform, sys
import importlib.metadata as im
import pandas as pd

OUT = Path(os.environ["INSPECTION_OUTPUT"])
OUT.mkdir(parents=True, exist_ok=True)
a = pd.read_csv("data/GSE242230_annotations.txt.gz", sep="\t")
assert not a.cell_id.duplicated().any()
a.groupby(["dataset", "tissue_type", "cell_type", "cell_type_specific", "filtered"], dropna=False).size().rename("n_cells").to_csv(OUT / "annotation_label_counts.tsv", sep="\t")
a.groupby(["sample_id", "cell_type_specific", "filtered"], dropna=False).size().rename("n_cells").to_csv(OUT / "sample_all_labels.tsv", sep="\t")
labels = ["Malignant - Basal", "Malignant - Classical", "Normal Epithelial"]
e = a.loc[(~a.filtered) & a.cell_type_specific.isin(labels)].copy()
c = pd.crosstab(e.sample_id, e.cell_type_specific).reindex(columns=labels, fill_value=0)
c["Malignant_total"] = c["Malignant - Basal"] + c["Malignant - Classical"]
c.to_csv(OUT / "epithelial_counts_by_sample.tsv", sep="\t")

# Read raw matrix headers only: no expression testing or expression QC yet.
rows = []
for p in sorted(Path("data").glob("*_matrix.mtx.gz")):
    prefix = str(p).removesuffix("_matrix.mtx.gz")
    sample = p.name.split("_")[1]
    with gzip.open(p, "rt") as f:
        header = next(f).strip()
        line = next(f)
        while line.startswith("%"):
            line = next(f)
        ng, nc, nnz = map(int, line.split())
    with gzip.open(prefix + "_features.tsv.gz", "rt") as f:
        features = f.read()
    fs = [line.split("\t") for line in features.splitlines()]
    expected = set(a.loc[a.sample_id.eq(sample), "cell_id"])
    seen = set()
    nbarcodes = 0
    with gzip.open(prefix + "_barcodes.tsv.gz", "rt") as f:
        for line in f:
            nbarcodes += 1
            key = sample + "_" + line.strip()
            if key in expected:
                seen.add(key)
    row = dict(sample_id=sample, matrix_format=header, n_features=ng,
               n_matrix_barcodes=nc, n_matrix_entries=nnz,
               n_feature_lines=len(fs), n_barcode_lines=nbarcodes,
               n_unique_feature_ids=len(set(r[0] for r in fs)),
               n_unique_symbols=len(set(r[1] for r in fs)),
               feature_types=";".join(sorted(set(r[2] for r in fs))),
               feature_sha256=hashlib.sha256(features.encode()).hexdigest(),
               n_annotated=len(expected), n_annotation_matches=len(seen),
               missing_annotation_ids=";".join(sorted(expected-seen)))
    rows.append(row)
m = pd.DataFrame(rows)
m.to_csv(OUT / "matrix_inventory_and_matches.tsv", sep="\t", index=False)

# Parse only supplied metadata; do not infer missing treatment information.
soft = Path("data/sample_metadata/GSE242230_sample_records.soft").read_text()
records = []
for block in soft.split("^SAMPLE = ")[1:]:
    lines = block.splitlines()
    r = {"geo_accession": lines[0].strip()}
    for line in lines[1:]:
        if line.startswith("!Sample_title = "):
            r["sample_id"] = line.split(" = ", 1)[1].split(",")[0]
        if line.startswith("!Sample_characteristics_ch1 = "):
            k, v = line.split(" = ", 1)[1].split(": ", 1)
            r[k] = v
    records.append(r)
s = pd.DataFrame(records)
s.to_csv(OUT / "soft_sample_metadata.tsv", sep="\t", index=False)
section = None
records = []
for line in Path("data/sample_metadata/Table_S1_sample_metadata.tsv").read_text().splitlines():
    fields = line.split("\t")
    if len(fields) < 5:
        continue
    patient = fields[1].strip()
    if patient in ["Surgically resected tumors", "EUS-FNB sample obtained"]:
        section = patient
    elif patient.startswith("WU") and patient[2:].isdigit():
        records.append(dict(sample_id="wu" + patient[2:].zfill(2), table_patient_id=patient,
                            section=section, table_stage=fields[2].strip(),
                            table_location=fields[3].strip(), table_pathology=fields[4].strip()))
t = pd.DataFrame(records)
s.merge(t, on="sample_id", how="outer", validate="one_to_one").to_csv(OUT / "joined_sample_metadata.tsv", sep="\t", index=False)
summary = {
    "n_annotation_rows": len(a), "n_retained_annotations": int((~a.filtered).sum()),
    "n_filtered_annotations": int(a.filtered.sum()),
    "n_epithelial_selected": len(e), "n_samples": len(m),
    "epithelial_counts": e.cell_type_specific.value_counts().to_dict(),
    "paired_samples_at_minimum_cells_per_group": {str(n): c.index[(c.Malignant_total >= n) & (c["Normal Epithelial"] >= n)].tolist() for n in [1, 5, 10, 20, 30]},
    "all_annotations_match_raw_barcodes": bool((m.n_annotated == m.n_annotation_matches).all()),
    "unique_feature_orders": int(m.feature_sha256.nunique()),
    "metadata_columns": s.columns.tolist(),
    "python": sys.version, "platform": platform.platform(),
    "packages": {d.metadata["Name"]: d.version for d in im.distributions()},
    "scope": "Annotation/metadata summaries, feature/barcode integrity, and matrix headers only. No DE, no matrix-entry loading."
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
print("\nEPITHELIAL COUNTS\n" + c.to_string())
print("\nMETADATA\n" + s.to_string(index=False))
print("\nMATRIX INVENTORY\n" + m.to_string(index=False))

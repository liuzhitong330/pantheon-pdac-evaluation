# Public data and downloads

Source: [GEO GSE242230](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242230) and the [associated study](https://www.nature.com/articles/s41698-023-00455-z). The repository includes derived results from this public dataset; it does not distribute the large original matrices.

For **comparison of saved results**, no GEO download is needed. Follow [REPRODUCING.md](REPRODUCING.md) to unpack the twenty saved comparison inputs.

For a **future raw-data replay**, download the series raw archive and original annotations into a new directory, from the repository root:

```sh
mkdir -p work/raw-data
curl --fail --location --retry 3 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242230/suppl/GSE242230_RAW.tar \
  --output work/raw-data/GSE242230_RAW.tar
curl --fail --location --retry 3 \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242230/suppl/GSE242230_annotations.txt.gz \
  --output work/raw-data/GSE242230_annotations.txt.gz
tar -xf work/raw-data/GSE242230_RAW.tar -C work/raw-data
mkdir -p work/raw-data/sample_metadata
cp data/metadata/GSE242230_sample_records.soft \
   data/metadata/Table_S1_sample_metadata.tsv work/raw-data/sample_metadata/
(cd work/raw-data && shasum -a 256 -c ../../data/SHA256SUMS.raw)
```

Linux users may substitute `sha256sum -c` for `shasum -a 256 -c`. Do not rename the 75 extracted files: the scripts identify the 25 samples from `GSM…_wu…_{matrix.mtx,features.tsv,barcodes.tsv}.gz`. The annotation file provides sample-prefixed cell IDs. The checksum list comes from the completed read-only input audit; a mismatch should be investigated before replay rather than silently accepted. Downloads and analysis are not run by cloning this repository.

Two small public metadata extracts are included unchanged:

| File | Origin and extraction |
|---|---|
| [GSE242230_sample_records.soft](../data/metadata/GSE242230_sample_records.soft) | All 25 sample blocks from the [GEO family SOFT file](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242230/soft/GSE242230_family.soft.gz), beginning at the first `^SAMPLE` record. |
| [Table_S1_sample_metadata.tsv](../data/metadata/Table_S1_sample_metadata.tsv) | Cells A1:N40 of sheet “Table S1” from the study's [supplementary workbook](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41698-023-00455-z/MediaObjects/41698_2023_455_MOESM2_ESM.xlsx), preserving blanks and row order. |

Their source URLs and downloaded-source hashes are recorded in [metadata_provenance.json](../data/metadata/metadata_provenance.json). Sample IDs, pathology and clinical fields originate in these public records; no private clinical data were added. Respect the original data and article terms and cite the study when reusing them.

“Normal Epithelial” is an author annotation within tumor-biopsy samples, not a healthy-donor cohort. Exact specific labels, Ensembl IDs and sample/patient mapping must be preserved. A missing gene test is not a zero effect.

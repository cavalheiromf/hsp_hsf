# Salmon RNA-seq Quantification Design

## 1. Objective

Quantify the 58 paired-end RNA-seq samples for *Setaria viridis*, *Sorghum bicolor*, *Triticum aestivum*, and *Glycine max* with Salmon, preserving transcript-level estimates and producing gene-level matrices suitable for later differential-expression analysis.

This stage ends with validated abundance matrices. Statistical normalization, exploratory analysis, differential-expression testing, and biological contrasts are separate downstream work.

## 2. Inputs

The workflow uses:

- 58 paired-end RNA-seq samples in `data/fastq/`, represented by 116 compressed R1/R2 FASTQ files;
- the curated sample metadata in `reference/external/SRA list_WD_eCO2 NCBI - Página1.csv`;
- the four Salmon 1.11.4 transcriptome indices in `results/rnaseq/<species>/index/`;
- the transcript FASTA files used to build those indices;
- the matching Ensembl Plants release 63 GFF3 annotations;
- a complete transcript-to-gene mapping derived from those GFF3 files.

Every sample must be associated with exactly one species, one R1 file, one R2 file, and the matching species index before it can be submitted.

## 3. Sample manifest

Create `config/rnaseq_samples.tsv` as the canonical machine-readable manifest. It contains at least:

```text
sample_id
species
bioproject
tissue
condition
replicate
layout
r1
r2
salmon_index
canary
```

The manifest is derived from the curated CSV rather than maintained independently. Validation fails on duplicate sample identifiers, unknown species, non-paired layouts, missing FASTQ files, or incomplete R1/R2 pairs.

## 4. Canary stage

The first quantification run contains one sample per species:

| Species | Canary | Rationale |
|---|---|---|
| *Triticum aestivum* | SRR33083344 | Original canary used to validate download and conversion |
| *Glycine max* | SRR26553587 | Small, well-watered sample |
| *Setaria viridis* | SRR39669459 | Ambient-CO2, well-watered control |
| *Sorghum bicolor* | SRR14935704 | Smallest available paired sample |

`SRR33083344` was processed before the main SRA array as the original canary. Its absence from the array log is expected and is not missing provenance; the manifest records its canary role explicitly.

The remaining 54 samples are not released until the four canary results have been reviewed.

## 5. Salmon quantification

A single Slurm array job supports both the canary subset and the complete manifest. Each task routes a sample to its species-specific index and runs mapping-based Salmon quantification on the compressed paired FASTQs.

The required Salmon settings are:

- automatic library detection with `-l A`;
- mapping validation with `--validateMappings`;
- sequence-bias correction with `--seqBias`;
- GC-bias correction with `--gcBias`;
- explicit thread count inherited from `SLURM_CPUS_PER_TASK`.

Trimming and inferential replicates are excluded from the initial run. Trimming is reconsidered only if canary diagnostics show adapter or read-quality problems. Decoy-aware indices are reconsidered if the existing transcriptome-only indices show poor or implausible mapping, particularly for hexaploid wheat.

Each sample writes to:

```text
results/rnaseq/quant/<species>/<sample_id>/
```

The output retains `quant.sf`, `meta_info.json`, `lib_format_counts.json`, auxiliary Salmon files, command parameters, and separate Slurm standard-output and standard-error logs.

## 6. Idempotence and failure handling

Before running Salmon, each task validates that:

1. R1 and R2 exist and are non-empty;
2. both gzip streams are readable;
3. the selected index exists and contains its completion metadata;
4. the output target does not contain an unreviewed partial run;
5. the filesystem has sufficient working space.

A completed sample is skipped only when its required Salmon outputs pass validation. Partial output is never silently overwritten. The job exits non-zero and reports the exact failing sample and check.

## 7. Transcript-to-gene mapping

Salmon quantifies every transcript present in the index. Gene-level summarization therefore requires a complete two-column `tx2gene` relation containing every indexed transcript and its parent gene.

The existing `results/rnaseq/gene_transcript_protein.tsv` is not used for this purpose because it contains only the representative protein-coding isoform selected for the HMMER pipeline. Using it would omit alternative transcripts.

One mapping is generated per species from the same Ensembl Plants release 63 GFF3 file used to build that species' transcript FASTA and Salmon index. Validation requires:

- unique transcript identifiers within each species;
- exactly one parent gene per mapped transcript;
- identifier compatibility between `tx2gene` and `quant.sf`;
- an explicit report of any indexed transcript without a gene mapping.

The mapping includes every indexed transcriptional feature with an `ID` and exactly one `Parent=gene:` relation, including both `mRNA` and non-coding RNA features such as tRNA, rRNA, and small RNAs.

Unmapped transcripts are not silently discarded. Any non-zero gap blocks gene-matrix publication until reviewed.

## 8. Matrix generation with tximport

The Bioconductor package `tximport` imports every Salmon `quant.sf`. Import and summarization are performed independently for each species; cross-species gene or transcript matrices are not created because their feature identifiers and experimental designs are not directly comparable. For each species, `tximport` produces:

- transcript-level estimated-count matrix;
- transcript-level TPM matrix;
- gene-level estimated-count matrix;
- gene-level TPM matrix;
- sample-specific effective-length information used by downstream statistical models;
- sample metadata aligned exactly to all matrix columns.

Gene-level values are summarized through the complete species-specific `tx2gene` mapping. The primary import uses `countsFromAbundance = "no"`, preserving Salmon's estimated counts and the sample-specific effective-length matrix used by `DESeqDataSetFromTximport`. TPM values are retained separately for abundance reporting and visualization. Later differential-expression analysis must use the `tximport` count and length information through an appropriate model such as DESeq2, not TPM values as raw counts.

Intermediate canary matrices are generated for validation. Final matrices are generated only when all intended samples have valid quantification outputs.

## 9. Quality-control report and release gate

The canary report records, per sample:

- input fragments;
- mapped fragments;
- mapping rate;
- inferred library type and compatibility;
- Salmon version and index identity;
- processing time;
- output-validation status;
- transcript-to-gene mapping coverage.

A mapping rate below 70% produces a warning. A rate below 50% prevents automatic release of the complete array until the result is investigated. These are operational triage thresholds, not biological acceptance criteria.

Release of the remaining 54 samples requires explicit review of all four canary records. Trimming, contamination investigation, or decoy-aware index reconstruction may be selected if the diagnostics justify it.

## 10. Tests

Automated tests cover:

- construction and schema validation of the sample manifest;
- sample-to-species and sample-to-index routing;
- missing, duplicate, or partial FASTQ pairs;
- complete GFF3 transcript-to-gene extraction;
- transcript identifier agreement between the index, `quant.sf`, and `tx2gene`;
- detection of partial or malformed Salmon output;
- deterministic ordering of samples in matrices;
- transcript- and gene-level `tximport` aggregation on a small fixture.

The four-canary run is the integration test for Slurm execution, real compressed FASTQs, the production indices, and final output collection.

## 11. Deliverables

```text
config/rnaseq_samples.tsv
jobs/salmon_quant.sbatch
scripts/build_rnaseq_sample_manifest.py
scripts/build_tx2gene.py
scripts/import_salmon_tximport.R
results/rnaseq/quant/<species>/<sample>/
results/rnaseq/qc/canary_quantification.tsv
results/rnaseq/qc/quantification.tsv
results/rnaseq/matrices/<species>/transcript_counts.tsv
results/rnaseq/matrices/<species>/transcript_tpm.tsv
results/rnaseq/matrices/<species>/gene_counts.tsv
results/rnaseq/matrices/<species>/gene_tpm.tsv
results/rnaseq/matrices/<species>/effective_lengths.tsv
results/rnaseq/matrices/<species>/sample_metadata.tsv
```

## 12. Out of scope

This design does not include:

- differential-expression models or contrasts;
- PCA, clustering, heatmaps, or batch correction;
- transcript-usage or splicing analysis;
- genome-aligned BAM production;
- promoter, motif, coexpression, or phylogenetic analyses.

Those analyses begin only after quantification matrices and experimental metadata have passed validation.

## 13. References

- Salmon documentation: <https://salmon.readthedocs.io/en/stable/salmon.html>
- Bioconductor tximport vignette: <https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html>

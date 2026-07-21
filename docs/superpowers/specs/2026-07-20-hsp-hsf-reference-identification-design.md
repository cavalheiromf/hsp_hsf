# HSP/HSF Reference and Identification Pipeline Design

## 1. Objective

Build a reproducible pipeline to identify six heat-shock protein families and the HSF transcription-factor family in plant proteomes, while preparing standardized genome references for downstream RNA-seq analysis.

The first implementation was a controlled pilot using *Setaria viridis* and *Sorghum bicolor*. These species already have HSF annotations in `reference/external/InterproScan.xlsx`, allowing the new pipeline to be validated against known results. The reference, HMMER, InterProScan, and catalog stages have since been expanded to *Triticum aestivum* and *Glycine max*.

## 2. Scope

### 2.1 Pilot scope

The completed pilot comprised:

1. download genome FASTA, GFF3, and complete protein FASTA files for *Setaria viridis* and *Sorghum bicolor*;
2. select one representative protein isoform per gene;
3. search the representative proteomes with seven Pfam HMM profiles;
4. extract candidate HSP/HSF proteins;
5. annotate candidates with the EMBL-EBI InterProScan REST service;
6. build per-species and combined catalogs;
7. compare newly detected HSFs with the existing HSF data.

### 2.2 Four-species expansion

After successful pilot validation, the same pipeline was extended to trigo (*Triticum aestivum*) and soja (*Glycine max*). Both species now have frozen references, representative proteomes, seven-family HMMER searches, InterProScan confirmation, and entries in the combined catalog. HSFs already present for the original 12 species are not replaced silently; new and discrepant classifications remain traceable.

Genome FASTA and GFF3 files are required only for the four RNA-seq species:

- *Triticum aestivum*;
- *Glycine max*;
- *Sorghum bicolor*;
- *Setaria viridis*.

## 3. Frozen reference source

All four RNA-seq references use Ensembl Plants release 63. A single reference assembly is used per species, even when multiple BioProjects exist for that species.

The selected assemblies are:

| Species | Assembly | INSDC accession |
|---|---|---|
| *Triticum aestivum* | IWGSC | GCA_900519105.1 |
| *Glycine max* | Glycine max v2.1 | GCA_000004515.4 |
| *Sorghum bicolor* | Sorghum bicolor NCBIv3 | GCA_000003195.3 |
| *Setaria viridis* | Setaria viridis v2.0 | GCA_005286985.1 |

The completed reference stage includes all four RNA-seq species. Every download records the source URL, Ensembl release, assembly, date, file size, and checksum.

## 4. Target protein families

The initial catalog is restricted to the seven families already defined in `reference/external/PFAM-info.txt`:

| Family | Pfam profile |
|---|---|
| Hsp20 | PF00011 |
| Hsp40/DnaJ | PF00226 |
| Hsp60 | PF00118 |
| Hsp70 | PF00012 |
| Hsp90 | PF00183 |
| Hsp100 | PF02861 |
| HSF | PF00447 |

Additional cochaperones and related families are outside the initial scope.

## 5. Repository layout

```text
config/
├── species.tsv
└── pfam_families.tsv

data/
├── reference/
│   └── ensembl_plants_63/
│       ├── setaria_viridis/
│       │   ├── genome.fa.gz
│       │   ├── annotation.gff3.gz
│       │   ├── proteins_all.fa.gz
│       │   └── checksums.tsv
│       └── sorghum_bicolor/
│           ├── genome.fa.gz
│           ├── annotation.gff3.gz
│           ├── proteins_all.fa.gz
│           └── checksums.tsv
└── pfam/
    ├── models/
    └── metadata.tsv

work/
├── representative_proteomes/
├── isoform_mapping/
├── hmmer/
└── candidates/

results/
├── interproscan/
├── catalogs/
└── validation/

scripts/
├── download_references.py
├── select_representative_isoforms.py
├── run_hmmer.sh
├── extract_candidates.py
├── classify_domains.py
└── validate_known_hsfs.py

jobs/
└── *.sbatch
```

Large downloaded and generated sequence files must remain outside Git. Configuration, scripts, manifests, validation reports, and final tabular catalogs should be versioned.

## 6. Data flow

```text
Ensembl Plants release 63
    -> genome FASTA + GFF3 + complete protein FASTA
    -> gene/transcript/protein mapping
    -> representative protein selection
    -> representative proteome
    -> seven independent hmmsearch runs
    -> union and extraction of candidates
    -> InterProScan annotation
    -> HSP/HSF catalog
    -> comparison with existing HSFs
```

Each intermediate stage must have a defined input, output, log, and validation check. Outputs from a failed or incomplete upstream stage must not be consumed.

## 7. Representative isoform selection

All protein isoforms are preserved in the raw downloaded FASTA. Exactly one translated protein is selected per protein-coding gene for the search proteome.

The selection order is:

1. select the transcript marked as canonical by Ensembl when such a designation is available;
2. otherwise select the longest translated protein;
3. if multiple proteins have the same maximum length, select the protein with the lexicographically smallest stable identifier.

The selection must be deterministic and produce a mapping table containing:

- gene identifier;
- selected transcript identifier;
- selected protein identifier;
- protein length;
- selection rule (`canonical`, `longest`, or `length_tie_id`);
- all alternative transcript and protein identifiers.

The pipeline must verify that each gene appears at most once in the representative proteome and that all selected protein identifiers exist in the downloaded FASTA.

## 8. HMMER candidate detection

The cluster-provided HMMER 3.4 module will be used. Searches are executed independently for each species and each of the seven Pfam profiles.

Pfam gathering thresholds are preferred through `hmmsearch --cut_ga`. This avoids defining a single arbitrary e-value threshold across heterogeneous protein families. The Pfam release and checksums of the HMM profiles must be recorded.

For every accepted match, retain:

- species;
- gene, transcript, and protein identifiers;
- family and Pfam accession;
- full-sequence score and e-value;
- domain score and conditional e-value;
- alignment and envelope coordinates;
- HMM and target coverage;
- source `domtblout` file.

A protein may initially match more than one profile. Such candidates are retained until domain architecture classification.

## 9. InterProScan confirmation

InterProScan will run through the EMBL-EBI Job Dispatcher REST service on the union of HMMER candidates, not on the complete proteomes. This avoids the large standalone installation during the two-species pilot. Inputs are split into deterministic batches of 100 proteins and submitted sequentially. Each batch records its sequence checksum, job identifier, status, timestamps, and result path so an interrupted run can resume without duplicating completed submissions. A valid contact email is supplied at runtime through `INTERPRO_EMAIL` and is never committed.

The requested InterProScan applications should include, when supported by the installed release:

- Pfam;
- SMART;
- PANTHER;
- CATH-Gene3D;
- SUPERFAMILY;
- PRINTS;
- PROSITE patterns and profiles.

TSV is the required machine-readable output. Additional JSON or GFF3 output may be retained if it adds no substantial execution cost.

InterProScan disagreements must be preserved instead of silently filtered. Each candidate receives one of these statuses:

- `confirmed`: expected HMMER and InterProScan signatures agree;
- `hmm_only`: HMMER match without the expected InterProScan confirmation;
- `interpro_related`: a related InterPro signature is present, but the expected signature is absent;
- `multi_family`: signatures support more than one target family;
- `review`: unusual architecture or unresolved evidence.

## 10. Catalog schema

The catalog must contain at least:

```text
species
gene_id
transcript_id
protein_id
protein_length
isoform_selection
family
pfam_id
hmm_score
hmm_evalue
domain_score
domain_evalue
domain_start
domain_end
domain_coverage
interpro_id
interpro_description
classification_status
reference_release
assembly_accession
```

Rows represent protein-family assignments. A protein with evidence for multiple families may have multiple rows and must also be flagged for review.

## 11. Existing HSF validation

The pilot will compare its HSF results with the historical PF00447 entries in `InterproScan.xlsx`. The workbook currently contains 23 unique *S. viridis* entries and 24 unique *S. bicolor* entries (47 total); the validation report records these observed counts rather than assuming the older 24/25 summary.

Validation proceeds in this order:

1. exact protein identifier match;
2. normalized identifier match for known formatting differences;
3. mapping through gene/transcript relationships in the release-63 GFF3;
4. sequence-based mapping using HSF sequences from `Alinhamento_HSF.fas` when identifiers changed.

Every known HSF must be mapped or assigned an explicit discrepancy category:

- unchanged identifier;
- updated identifier;
- different representative isoform;
- absent from the release-63 proteome;
- not recovered by the HMM search;
- sequence or annotation inconsistency.

Potential new HSF candidates must be reported separately rather than merged silently into the historical set.

## 12. Slurm execution model

The workflow consists of dependent stages:

```text
01_download_references
    -> 02_select_isoforms
    -> 03_hmmsearch_array
    -> 04_extract_candidates
    -> 05_interproscan
    -> 06_build_catalog
    -> 07_validate_known_hsfs
```

Dependencies use Slurm `afterok`, ensuring that downstream stages do not run after an upstream failure. HMMER tasks should use a species-by-family job array. InterProScan uses a Slurm array limited to one active species (`%1`); within each species, API batches are submitted sequentially and polled every 30 seconds to respect the EMBL-EBI fair-use policy.

Each task must write separate standard-output and standard-error logs containing the species, family or stage, tool version, command line, start time, end time, node, and exit status.

## 13. Validation and failure handling

The pipeline must fail early when any of the following checks fails:

- missing or mismatched checksum;
- corrupt gzip file;
- empty FASTA or GFF3;
- missing GFF3-to-FASTA identifier mapping;
- empty representative proteome;
- duplicated selected genes;
- missing Pfam model;
- empty HMMER output for all families;
- missing InterProScan output;
- malformed catalog schema.

Partial or ambiguous biological results do not cause silent deletion. They are retained with a review status and summarized in the validation report.

## 14. Pilot acceptance criteria

The pilot is successful when:

1. release-63 reference files for both species are downloaded, verified, and documented;
2. the representative-proteome selection is deterministic and produces no duplicate genes;
3. all seven Pfam searches complete with traceable outputs;
4. InterProScan completes on the candidate sets;
5. preliminary catalogs for the six HSP families and HSF are produced;
6. all historical PF00447 entries present in the workbook are mapped or assigned an explicit discrepancy explanation;
7. commands, versions, parameters, logs, and checksums are sufficient to reproduce the run;
8. the pipeline can be expanded by adding species to configuration rather than rewriting scripts.

## 15. Explicitly deferred work

The following work is not part of the pilot:

- rerunning HSF identification for all 12 previously annotated species;
- identifying cochaperone families outside the seven specified Pfam models;
- statistical analysis of quantified RNA-seq reads, which is specified separately after the Salmon quantification stage;
- differential-expression analysis;
- phylogenetic tree reconstruction;
- promoter, motif, or coexpression analysis;
- genome-aligned BAM production for wheat, soybean, sorghum, or Setaria.

RNA-seq quantification is defined in `2026-07-21-salmon-rnaseq-quantification-design.md`. The other stages depend on successful completion and review of the four-species quantification outputs.

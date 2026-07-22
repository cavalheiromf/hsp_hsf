# RNA-seq Exploratory QC and Experimental Design

**Date:** 2026-07-22

## Objective

Create a reproducible post-quantification quality-control workflow for the four plant species in this project. The workflow will normalize experimental metadata, assess sample-level behavior, document estimable differential-expression contrasts, and publish an HTML report plus reusable TSV and image files. Differential-expression testing is explicitly outside this scope.

## Inputs

The workflow consumes only versioned, validated products of the Salmon quantification pipeline:

- `config/rnaseq_samples.tsv`
- `results/rnaseq/qc/quantification.tsv`
- `results/rnaseq/matrices/final/<species>/gene_counts.tsv`
- `results/rnaseq/matrices/final/<species>/sample_metadata.tsv`

The other final matrices remain available for audit but are not required for exploratory normalization. Sample identifiers and their order must agree exactly across the manifest, QC table, metadata, and count matrices.

## Scope and Boundaries

The workflow will:

- create explicit biological factors from the free-text `condition` field;
- preserve the original condition text for provenance;
- analyze each species independently;
- analyze wheat and soybean BioProjects independently while also providing a joint, batch-colored visualization for each species;
- normalize counts for exploratory visualization;
- produce library, detection, correlation, distance, and PCA diagnostics;
- retain every sample and record advisory flags with reasons;
- document statistically estimable contrasts without fitting differential-expression models;
- publish HTML, TSV, PNG, and SVG artifacts.

The workflow will not:

- remove or overwrite quantified Salmon results;
- exclude samples automatically;
- perform hypothesis testing for differential expression;
- compare raw expression values directly across species;
- merge independent BioProjects into a single biological contrast.

## Architecture

### Metadata normalization

`scripts/build_rnaseq_design.R` will read the canonical manifest and produce `config/rnaseq_design.tsv`. Each row will retain `sample_id`, `species`, `bioproject`, `tissue`, `condition`, `replicate`, and `canary`, and add normalized factors including `analysis_group`, `co2`, `water`, and `dehydration_stage` where applicable.

The parser will use explicit mappings for the current BioProjects rather than general keyword guessing. Any unrecognized condition will be a blocking error. Factor reference levels will be recorded in the output and contrast table.

### Exploratory calculations

`scripts/run_rnaseq_exploratory_qc.R` will validate inputs, build one exploratory object per analysis group, calculate metrics, and generate reusable tables and figures. It will use Bioconductor DESeq2 for size-factor estimation and variance-stabilizing transformation and tidyverse/ggplot2 for tabulation and plots.

### Reproducible report

`reports/rnaseq_exploratory_qc.qmd` will consume the generated tables and figures and render `reports/rnaseq_exploratory_qc.html`. The report will describe methods, summarize every analysis group, show advisory flags, and list possible contrasts. It will not contain hidden analytical state required to reproduce the results; the scripts remain the computational interface.

## Experimental Designs

| Species | BioProject | Analysis group | Biological factors | Replicates | Estimable model |
|---|---|---|---|---:|---|
| *Setaria viridis* | PRJNA1496374 | `setaria_prjna1496374` | ambient/elevated CO2; control/cycle 1/cycle 3 dehydration | 3 | `~ co2 * dehydration_stage` |
| *Sorghum bicolor* | PRJNA742236 | `sorghum_prjna742236` | ambient/elevated CO2 under drought | 4 | `~ co2` |
| *Triticum aestivum* | PRJNA1249400 | `wheat_prjna1249400` | ambient/elevated CO2; control/drought | 3 | `~ co2 * water` |
| *Triticum aestivum* | PRJNA793265 | `wheat_prjna793265` | control/drought under elevated CO2 | 3 | `~ water` |
| *Glycine max* | PRJNA295411 | `soybean_prjna295411` | ambient/elevated CO2; well-watered/air-drying | 2 | `~ co2 * water` |
| *Glycine max* | PRJNA1033144 | `soybean_prjna1033144` | control/drought under elevated CO2 | 3 | `~ water` |

`dehydration_stage` is categorical. Cycle 1 and cycle 3 will not be treated as numeric time points. For the two factorial designs, the contrast inventory will include estimable simple effects and interaction terms with explicit numerator, denominator, and reference levels.

## Statistical Processing

### Count filtering

Filtering is performed separately within each analysis group. A gene is retained when its estimated count is at least 10 in at least as many samples as the smallest biological replicate group in that analysis group. The filter status and retained-gene totals will be written to tables.

### Exploratory normalization

Salmon gene counts may be fractional. They will be rounded only when constructing the DESeq2 object. Size factors will be estimated within each analysis group. A blind variance-stabilizing transformation (`blind = TRUE`) will be used for unsupervised PCA, sample correlation, and distance diagnostics. These transformed values are for QC only and will not be presented as differential-expression results.

### Metrics and visualizations

For every analysis group, the workflow will produce:

- raw library size and normalized size factor;
- number of genes detected at counts of at least 1 and at least 10;
- number and proportion of genes retained after filtering;
- mapping rate and Salmon library-format information from the consolidated QC table;
- VST expression distributions;
- PCA scores and variance explained;
- sample-to-sample Spearman correlation;
- Euclidean sample distance on VST values;
- heatmaps of correlation and distance;
- plots colored by biological factors and labeled by sample.

For wheat and soybean, an additional species-level PCA will be produced solely to visualize BioProject separation. It will be colored by BioProject and must not be used to define a cross-project biological contrast.

### Advisory flags

No sample is removed. Flags will be emitted as a table containing `sample_id`, metric, observed value, rule, severity, and explanation. The consolidated Salmon mapping warnings for `SRR39669466` and `SRR39669467` will be preserved. Additional distributional diagnostics may flag samples for review, but the report must distinguish an advisory flag from a blocking input-validation error.

Because several biological groups contain only two or three replicates, PCA position alone will never automatically define an outlier. The report will show ranks and diagnostic values so that exclusions, if ever considered later, require an explicit biological decision and a new analysis.

## Contrast Inventory

`results/rnaseq/exploratory_qc/contrasts.tsv` will document, without executing, the possible DESeq2 analyses. Required columns are:

- `contrast_id`
- `species`
- `bioproject`
- `analysis_group`
- `model_formula`
- `effect_type`
- `numerator`
- `denominator`
- `reference_levels`
- `estimable`
- `notes`

Simple one-factor projects will contain the direct treatment contrast. Factorial projects will contain main effects at stated reference levels, simple effects, and the interaction term. Every row must be derivable from observed design cells and biological replicates.

## Outputs

The workflow will publish:

```text
config/rnaseq_design.tsv
results/rnaseq/exploratory_qc/
  contrasts.tsv
  session_info.txt
  tables/
    sample_metrics.tsv
    gene_filter_summary.tsv
    pca_scores.tsv
    advisory_flags.tsv
  figures/
    <analysis_group>_library_size.{png,svg}
    <analysis_group>_detected_genes.{png,svg}
    <analysis_group>_vst_distribution.{png,svg}
    <analysis_group>_pca.{png,svg}
    <analysis_group>_correlation_heatmap.{png,svg}
    <analysis_group>_distance_heatmap.{png,svg}
    wheat_joint_bioproject_pca.{png,svg}
    soybean_joint_bioproject_pca.{png,svg}
reports/rnaseq_exploratory_qc.qmd
reports/rnaseq_exploratory_qc.html
```

Additional per-group matrix-like TSV files may be stored under `tables/` when needed to reproduce a heatmap, but normalized expression matrices will not be published unless a downstream interface requires them.

## Validation and Failure Handling

The workflow blocks publication when it encounters:

- missing or duplicate samples;
- matrix columns that differ from manifest order;
- unrecognized conditions or factor combinations;
- absent or duplicated replicate identifiers within a design cell;
- nonfinite or negative counts;
- incompatible matrix dimensions;
- a biological cell without the documented number of replicates;
- unavailable required R or Quarto dependencies.

All report products will be assembled in a sibling staging directory. The staged collection must pass schema, dimension, finiteness, file-inventory, and metadata-alignment checks before promotion. If processing or validation fails, the previously published collection and HTML report remain unchanged. Promotion will use same-filesystem renames with rollback of the prior collection.

## Testing

Tests will use small synthetic count matrices and manifests. They must cover:

- exact factor parsing for every current BioProject;
- rejection of unknown condition strings and invalid replicate cells;
- correct reference levels and contrast inventory;
- group-specific low-count filtering;
- exact sample and feature alignment;
- rejection of negative and nonfinite counts;
- preservation of all samples carrying advisory flags;
- deterministic table schemas and sample order;
- expected figure inventory;
- successful Quarto rendering;
- late-failure preservation of an existing published report;
- cleanup of staging and backup artifacts after success.

The production acceptance check requires all expected groups, all 58 samples in metrics, no silent exclusions, both accepted Setaria mapping warnings, and an HTML report that renders without errors.

## Reproducibility

`session_info.txt` will record R, Bioconductor, DESeq2, tidyverse, ggplot2, heatmap package, and Quarto versions. All input and output paths in scripts and the Quarto source will be relative to the repository root. Generated tables will use UTF-8, tab delimiters, stable ordering, and explicit missing-value handling.

## Success Criteria

The work is complete when:

1. `rnaseq_design.tsv` contains all 58 samples with validated factors and analysis groups.
2. Six documented analysis groups pass design validation.
3. All samples appear in the QC metrics and report; none are automatically removed.
4. The accepted Setaria mapping warnings remain visible.
5. Every expected TSV, PNG, SVG, and HTML file is present and validated.
6. The contrast inventory contains only observed and estimable comparisons.
7. Synthetic tests and the production render pass.
8. A failed rerun cannot replace a previously valid report collection.

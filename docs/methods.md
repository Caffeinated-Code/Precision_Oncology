# Methods

## Cohort

The primary analysis uses the public cBioPortal study `skcm_dfci_2015`, corresponding to Van Allen et al., Science 2015. The cohort contains whole-exome sequencing from metastatic melanoma tumor-normal pairs treated with CTLA-4 blockade.

## Response Grouping

The cBioPortal clinical field `DURABLE_CLINICAL_BENEFIT` is converted into a binary analysis label:

- durable clinical benefit: `CR`, `PR`, `SD`
- no durable benefit: `PD`
- excluded from binary response analysis: `X`, missing, or indeterminate values

This grouping follows the common immuno-oncology convention of treating complete response, partial response, and stable disease as benefit categories, while preserving uncertainty for indeterminate labels.

## Mutation Processing

The analysis downloads mutation records from the cBioPortal mutation profile `skcm_dfci_2015_mutations`. Silent and noncoding annotation categories are excluded from the main screen. A gene is counted as mutated in a sample if at least one nonsynonymous mutation is observed for that gene in that sample.

The primary unit for gene enrichment is therefore sample-level mutation presence, not raw variant count.

## Gene Enrichment

Each gene is tested using a 2 x 2 Fisher exact test:

| | Mutated | Wild-type |
|---|---:|---:|
| Durable clinical benefit | a | b |
| No durable benefit | c | d |

Benjamini-Hochberg q-values are reported across tested genes.

Frequency tiers:

- `primary_screen`: gene mutated in at least 5 analyzed samples
- `low_frequency`: gene mutated in 2-4 analyzed samples
- `singleton`: gene mutated in 1 analyzed sample and not tested in the main table

Low-frequency findings are reported separately because sparse events can create unstable odds ratios.

## Burden And Immune Context

The script summarizes available cBioPortal clinical fields for mutation count, nonsynonymous mutation load, TMB, neoantigen load, and selected CIBERSORT/ESTIMATE immune features. Group differences are summarized using medians and rank-biserial effect sizes.

## Caveats

- The analysis is retrospective and exploratory.
- It does not claim clinical actionability.
- It does not distinguish clonal from subclonal mutations.
- It does not perform external validation.
- Gene-level mutation presence does not imply functional impact.
- cBioPortal-derived values may not exactly match every preprocessing choice in the original publication.


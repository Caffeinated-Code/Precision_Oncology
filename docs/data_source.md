# Data Source

## Primary Public Cohort

- Study: Metastatic Melanoma (DFCI, Science 2015)
- cBioPortal ID: `skcm_dfci_2015`
- Citation: Van Allen et al., Science 2015
- PMID: 26359337
- Data portal: https://www.cbioportal.org/study/summary?id=skcm_dfci_2015
- Reference genome listed by cBioPortal: hg19

## Data Pulled By The Pipeline

The analysis script downloads:

- study metadata
- patient-level clinical annotations
- sample-level clinical annotations
- mutation calls from the study mutation profile

The generated files under `data/derived/` are derived from public cBioPortal API responses.

## Why This Cohort

This cohort is a strong precision-oncology case study because it connects:

- whole-exome somatic mutation profiling
- immune checkpoint therapy response
- tumor mutation and neoantigen burden
- candidate gene response associations
- immune microenvironment estimates available through cBioPortal

That combination supports a realistic translational analysis using public clinical-genomic data.

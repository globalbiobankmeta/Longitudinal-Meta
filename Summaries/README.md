# Cohort Summaries

This directory contains scripts to compute biobank-specific cohort summaries, such as sample size, follow up time, and age distributions. This relies on existing phenotype files used for GWAS, as well as basic information on enrollment and last followup. Runtime is expected to be a few seconds per phenotype. 

``gbmi_cohort_summaries.R`` is an all-in-one script to summarize overall biobank cohort, onset phenotypes, and progression phenotypes. Examples for running each analysis and sharing results are detailed below

- [PDF documentation](gbmi_cohort_summaries.pdf)
- [Download HTML documentation](https://raw.githubusercontent.com/globalbiobankmeta/Longitudinal-Meta/main/Summaries/gbmi_cohort_summaries.html)

Latest updates to ``gbmi_cohort_summaries.R``

- Aug 28, 2026: modified to remove NAs by default from summaries (i.e., mean/median/sd); NAs can remain in the histogram, and otherwise can be kept by adding the optional flag ``--na-keep`` to the R command
- Aug 19, 2026: added ``--km-round`` flag as an option to round event times to nearest integer

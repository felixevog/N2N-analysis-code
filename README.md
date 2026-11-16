# Narrative-to-Network (N2N): Analysis Code

This repository contains R scripts associated with the manuscript:

**From Brief Narratives to Person-Specific Symptom Networks: Formalizing Subjective Models of Psychopathology**

## Overview

The study investigates whether brief prompted narratives can be coded into person-specific directed symptom networks, agreement between independent coders, correspondence with systematically elicited Perceived Causal Networks (PECAN), and exploratory associations with psychological distress.

The research includes a development phase (n = 22) and an independent evaluation phase (n = 50).

## Scripts

| File | Purpose |
| --- | --- |
| `N2N_analysis_pilotstudy.R` | Development-phase analyses and exploratory outputs. |
| `N2N_analysis_mainstudy.R` | Evaluation-phase analyses, including comparisons between narrative-derived and PECAN networks. |
| `N2N_analysis_combined_DASS.R` | Exploratory pooled analyses of network metrics and DASS-21 distress scores across phases. |
| `loop_visualisation_pilotstudy.R` | Visualizations of development-phase narrative and PECAN networks. |
| `loop_visualisation_mainstudy.R` | Visualizations of evaluation-phase narrative and PECAN networks. |

## Software

Analyses were conducted in **R 4.5.0**. The scripts reference packages including `tidyverse`, `psych`, `qgraph`, `PECAN2`, `irr`, `afex`, `emmeans`, `ppcor`, `visNetwork`, `htmlwidgets`, and `webshot2`. Dependencies vary by script. Package versions have not been locked in a reproducible environment.

## Data access and reproducibility

**The repository contains analysis scripts, not participant data.** Raw narratives and participant-level survey or network data are not uploaded because they may contain sensitive information and are subject to ethical and data protection restrictions. De-identified derived data may be available from the corresponding author upon reasonable request, subject to applicable restrictions.

The scripts expect inputs in a local project structure (notably `data/rawdata/` and `data/bearbeitet/`) and write outputs to `output/`. They **cannot be run end-to-end from this repository alone** without the relevant input files. The original scripts also contain ancillary, exploratory, and development-phase analyses beyond the subset presented in the manuscript.

The corresponding author's contact information is provided in the manuscript.

## Interpretation

The scripts are shared to improve transparency. Some sections reflect earlier methodological choices (including interrater agreement statistics), while the manuscript describes the final reported analyses, including Dice-based agreement. The availability of these files should **not** be interpreted as a verified, self-contained reproduction package for every published table or figure.

## License

No software license has been selected yet. Until a license is added, standard copyright restrictions apply.

## Threshold-Dependent Gene Drive Simulations in *Plutella xylostella*

This repository contains all code, data structures, and visualisation tools used in a high-resolution simulation study evaluating the suppression potential of a MEREA-based threshold-dependent gene drive system in *Plutella xylostella* (diamondback moth), a globally pervasive lepidopteran pest.

### Project Overview

The diamondback moth causes over $4–5 billion USD in annual agricultural losses due to its resistance to conventional insecticides. This project assesses the viability of **MEREA** (Maternal Effect Recessive Embryonic Arrest) gene drive systems, both one-locus and two-locus (split-drive), for population suppression under ecologically realistic scenarios.

Simulations were performed using **MGDrivE** in R, incorporating overlapping generations, persistent wild-type larval reservoirs, and custom inheritance matrices reflecting toxin-antidote logic. Results demonstrate that success is highly conditional, dependent on threshold frequencies, staggered releases, and absence of resistance alleles.

### Key Research Questions

- How does ecological realism (e.g. pre-existing larvae) affect invasion thresholds?
- What release regimes achieve successful suppression?
- How do one-locus and two-locus MEREA systems compare?
- What is the system's vulnerability to resistance emergence?

### Tools & Technologies

- **Language:** R (v4.4.2)
- **Core Packages:** [MGDrivE](https://cran.r-project.org/web/packages/MGDrivE/index.html), `ggplot2`, `tidyverse`, `stringr`
- **Simulation Type:** Deterministic, non-spatial
- **Code Assistance:** GPT-4 (used for debugging and code scaffolding)

### Reproducibility

All simulation scripts are contained in the `/scripts` directory and are fully reproducible.

Custom inheritance cubes for the MEREA systems are located in `/cubes`, with filenames indicating drive type (e.g. `Cube_MEREA_with_resistance_allele.R`, `Cube_MEREA_two_loci.R`).

### Key Findings

| Configuration | Threshold | Releases | Resistance | Suppression Achieved | Time to Collapse |
|---------------|-----------|----------|------------|----------------------|------------------|
| One-locus     | ≥50%      | ≥3       | Absent     | Yes                  | ~900 days        |
| Two-locus     | ≥90%      | 8        | Absent     | Yes                  | ~450 days        |
| Any           | Any       | Any      | 0.01%      | No                   | N/A              |

- **Single releases fail** under ecological realism due to larval carryover.
- **Staggered releases** are critical to surpassing invasion thresholds.
- **Resistance alleles** at just 0.01% fully blocked suppression.

### Repository Structure

```
📁 scripts/               # Main R scripts for simulation and visualisation
📁 cubes/                 # MEREA gene drive inheritance matrices
📁 mgdrive/               # Simulation data
📁 plots/                 # Generated visualisations for the report
📁 renv/                  # renv environment metadata
📄 renv.lock              # Exact package versions for reproducibility
📄 .Rprofile              # Initialises renv on project load
📄 .gitignore             # Files and folders excluded from version control
📄 README.md              # Project overview
```

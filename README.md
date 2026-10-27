# EvoTempo

**EvoTempo** is an interactive R Shiny toolkit for evaluating temporal signal in measurably evolving populations. It connects exploratory diagnostics with preparation and interpretation of external BEAST analyses.

EvoTempo provides these analyses within a unified graphical interface, reducing the need for complex command-line workflows while supporting reproducible analysis and publication-ready visualization.

## Modules and data connections
---
| Module | Function | Input and execution |
| --- | --- | --- |
| 1 Data Processing and Mantel Test | Parse dates and examine temporal–genetic association | One tree OR nucleotide FASTA alignment, with dates from identifiers or CSV |
| 2 Interactive Root-to-Tip Regression | Explore global, clade and discrete-trait regressions | Tree and dates shared from Module 1 |
| 3 Date-Randomization Test (DRT) | Generate randomized XMLs and compare rate posteriors | Supported BEAST XML or original/randomized MCMC logs; independent entry point |
| 4 Bayesian Evaluation of Temporal Signal (BETS) | Compare supplied log marginal likelihoods | External estimates for heterochronous and isochronous models; independent entry point | 

Module 1 supplies Module 2 after **Run Mantel Test** or **Send to RTT (Skip Mantel Test)**. Module 2 has no independent tree/date upload. Cluster assignments from Module 1 must be exported and uploaded explicitly into Module 3. The Mantel result does not change the **Enable Cluster-based DRT** setting.

BEAST execution, MCMC convergence assessment and marginal-likelihood estimation take place outside EvoTempo.

## Key Features & Modules

EvoTempo consists of four complementary analytical modules.

### 1. Data Parsing & Mantel Test

Evaluates potential **temporal–genetic confounding** that may generate an apparent temporal signal.

Key functions include:

* Sampling-date parsing and formatting.
* Neighbor-Joining (NJ) rooting and genetic-distance calculation.
* Permutation testing based on Pearson's correlation coefficient.
* Detection of significant temporal–genetic associations.
* Generation of **Cluster Metadata** for cluster-based DRT when confounding is detected.

### 2. Root-to-Tip Regression

Provides interactive exploratory assessment of temporal structure in phylogenetic trees.

Key functions include:

* Root optimization using NVE-RMS or maximum correlation (`r`).
* Global root-to-tip regression under a strict-clock framework.
* Interactive clade selection for local-clock analyses.
* Discrete-trait partitioning.
* Permutation-based assessment of the regression **R²** statistic.
* Export of vector plots and RTT plot data.

### 3. Date-Randomization Test (DRT) Framework

Provides an integrated workflow for generating and evaluating date-randomization analyses.

The module supports:

* **Standard permutation**
* **Cluster-based permutation**
* BEAST 1.x and BEAST 2.x XML configuration generation.
* MCMC log analysis and molecular-clock rate extraction.
* Burn-in exclusion and 95% HPD calculation.
* Evaluation according to the **CR1** and **CR2** criteria.

### 4. Bayesian Evaluation of Temporal Signal (BETS)

Evaluates temporal signal through Bayesian comparison of heterochronous and isochronous models.

EvoTempo supports:

* Comparison of **M<sub>het</sub>** and **M<sub>iso</sub>** models.
* Extraction of marginal likelihood estimates from BEAST analyses.
* Calculation of **log Bayes Factors**.
* Interpretation of Bayes Factors according to the framework of Kass & Raftery (1995).

---

## Workflow

EvoTempo follows a modular workflow in which the analytical modules can be used independently or connected through intermediate outputs.

```text
                    Sequence / Tree
                  + Sampling Metadata
                           │
          ┌────────────────┼────────────────┐
          │                │                │
          ▼                ▼                ▼
     Module 1          Module 2          Module 3
     Mantel Test       Root-to-Tip           DRT
          │             Regression            │
          │                                    │
      P < 0.05?                    ┌───────────┴───────────┐
          │                        ▼                       ▼
          ▼                   Standard DRT          Cluster-based DRT
   Cluster Metadata                                    ┌────┴────┐
          │                                            ▼         ▼
          └───────────────────────────────────────►  User CSV  Mclust
                                                     (from M1)  auto-cluster


                 Independent BEAST Analysis
                           │
                           ▼
                 Marginal Likelihoods
                           │
                           ▼
                      Module 4
                        BETS
```

The workflow supports both independent module-level analyses and integrated temporal-signal assessment, with **Cluster Metadata** enabling cluster-based date randomization when significant temporal–genetic confounding is detected.

---

## Project Structure
```text
EvoTempo/
├── app.R                       # Main application entry point
├── global.R                    # Global setup, package loading, and shared functions
├── R/                          # Backend algorithms and Shiny modules
│   ├── mod01_Mantel.R          # Tab 1: Data parsing, Mantel test, and cluster visualization
│   ├── mod02_RTT.R             # Tab 2: Interactive tree and root-to-tip regression
│   ├── mod03_DRT.R             # Tab 3: DRT workflow and MCMC log evaluation
│   ├── mod04_BETS.R            # Tab 4: Bayesian Evaluation of Temporal Signal
│   ├── mod05_About.R           # Tab 5: About, citations, and terms of use
│   ├── utils_xml.R             # BEAST XML parsing and date randomization
│   ├── utils_date_parsing.R    # Sampling-date extraction and decimal-year conversion
│   ├── utils_temporal.R        # Phylogenetic distances and root-to-tip analysis
│   └── utils_coda.R            # MCMC log processing and HPD/CR1/CR2 evaluation
├── inst/
│   └── extdata/                # Example datasets for testing
├── www/                        # Static web resources
│   └── logo.png                # EvoTempo logo
├── DESCRIPTION                 # Project metadata and R package dependencies
├── LICENSE                     # MIT license metadata
└── README.md                   # Project overview and usage guide
```
## Installation and launch

The declared minimum is R 4.1.0. Install the listed packages in a compatible R environment. `ggtree` comes from Bioconductor; `mclust` is needed only for automatic date clustering.

```r
install.packages(c(
  "ade4", "ape", "bslib", "coda", "crosstalk", "dplyr", "DT",
  "future", "ggplot2", "lubridate", "plotly", "promises", "seqinr",
  "shiny", "shinyjs", "stringr", "zip"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ggtree")
# Optional, for date-based automatic clustering:
install.packages("mclust")
setwd("path/to/EvoTempo") # directory containing app.R and global.R
shiny::runApp()
```

Specific minimum browser/OS versions have not been established for this snapshot. This directory contains a Shiny application; the `DESCRIPTION` file lists metadata and dependencies.

## Input & Output

### Input

Depending on the selected module, EvoTempo can work with:

* Sequence data and/or phylogenetic trees.
* Sampling-date metadata.
* Discrete trait metadata, such as geographic location or host.
* BEAST 1.x or BEAST 2.x XML configurations.
* BEAST MCMC log files.

### Output

EvoTempo can generate or export:

* Cluster Metadata for cluster-based DRT.
* BEAST XML configurations.
* Root-to-tip vector plots.
* RTT plot data.
* DRT results and 95% HPD intervals.
* Bayesian model-comparison results and Bayes Factors.

---

## External Software

EvoTempo is designed to work with **BEAST 1.x and BEAST 2.x** for Bayesian phylogenetic analyses required by the DRT and BETS workflows.

BEAST analyses are performed externally. EvoTempo generates the required XML configurations and subsequently processes the resulting MCMC logs for downstream temporal-signal evaluation.

---

## Citations

EvoTempo integrates established methodologies from molecular epidemiology, phylogenetics, and phylodynamics. If you use EvoTempo in your research, please cite the software together with the original publications relevant to the analyses performed.
Citation

If you use EvoTempo in your research, please cite:

**EvoTempo: An interactive R Shiny toolkit for temporal signal evaluation in measurably evolving populations.**

**Citation details**: Li et al. EvoTempo: An interactive R Shiny toolkit for temporal signal evaluation in measurably evolving populations. 2026; Manuscript in preparation.

Please also cite the relevant original methods and software used in your analysis.

---

## License

EvoTempo is distributed under the **MIT License**. See the [LICENSE](LICENSE) file for details.

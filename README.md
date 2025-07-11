# CRC Meta-Analysis II

## Meta-analysis reveals a universal microbiome signature for colorectal cancer irrespective of onset age and sequencing method

**Selin Pekel\*, Nicolai Karcher\*, Morgan Essex, Stefano Romano, Quinten R. Ducarmon, Fabian Springer, Christian Schudoma, Martin Larralde, Alberto Lupatin, Sebastian Zeissig, Michael Zimmermann, Georg Zeller#**\
(\*equal contribution, #correspondence: [georg.zeller\@gmail.com](mailto:georg.zeller@gmail.com){.email})

------------------------------------------------------------------------

## Overview

This repository contains the code, processed data, and figure scripts for the largest single-disease gut microbiome meta-analysis to date, integrating whole-genome shotgun (WGS) and 16S rRNA sequencing data from **27 studies (N = 6,779)** of colorectal cancer (CRC).

------------------------------------------------------------------------

## Requirements

-   This project requires **R version 4.4.1** or later.\
    \
    To ensure all required packages are installed, run the following command from the `crc-meta-ii/src/` directory:

    <div>

    Rscript requirements.R

    </div>

-   This project requires the SIAMCAT package, specifically version 2.7.2, as it includes important updates necessary for compatibility with the current workflow. Using an older version may lead to errors or unexpected behavior. To ensure proper functionality, please download the **SIAMCAT v2.7.2**

-   In addition, SHAP value calculations in this project were performed using classifiers built with the **mlr3 package (version 0.20.2)**, so please ensure that this version is installed to avoid compatibility issues.

------------------------------------------------------------------------

## Repository Structure

``` text
crc-meta-ii
├── data/
│   ├── Metadata.all.samples.tsv
│   ├── Metadata.all.samples.balanced.tsv
│   ├── Relab.all.samples.tsv
│   ├── Relab.all.samples.balanced.tsv
│   ├── Raw.counts.all.samples.tsv
│   ├── Raw.counts.wgs.motus.all.samples.tsv
│   ├── mOTUs_taxonomy_table.tsv
│   ├── fiber_data/
│   ├── functional_data/
│   └── results/
│       ├── scv.loso/
│       │   ├── ad.loso.test/
│       │   ├── crc.loso.test/
│       │   ├── crc.loso.train/
│       │   └── crc.scv.sst/
│       └── shap.analysis/
│           ├── fold_info/
│           ├── kernelshap_objects/
│           └── models/
├── figures/
│   ├── figure1/
│   ├── figure2/
│   ├── figure3/
│   ├── figure4/
│   ├── figure5/
│   ├── figure6/
│   ├── extended.data.figure1/
│   ├── extended.data.figure2/
│   ├── extended.data.figure3/
│   ├── extended.data.figure4/
│   ├── extended.data.figure5/
│   └── extended.data.figure7/
├── models/
│   ├── 01.Train.16s.wgs.rf.models.R
│   ├── 02.Train.eo.lo.crc.rf.models.R
│   ├── 03.Train.unified.crc.model.R
│   ├── 04.run.Train.LOSO.rf.models.sh
│   ├── 05.Train.SCV.SST.rf.models.R
│   ├── 06.Train.ad.ctr.model.R
│   ├── input.LOSO.rf.models.sh
│   └── Train.LOSO.rf.models.R
├── src/
│   ├── analysis/
│   │   ├── 01.PCoA.analysis.R
│   │   ├── 02.lmm.16S.WGS.R
│   │   ├── 03.lmm.EO.LO.R
│   │   ├── 04a.prepare.data.and.models.for.shap.analysis.R
│   │   ├── 04b.run_SHAP_test.sh
│   │   ├── 04c.run_SHAP_test.slurm
│   │   ├── 04d.get.shap.median.values.R
│   │   ├── 05.alpha.diversity.and.sequencing.depth.R
│   │   ├── 06.lmm.functional.profiles.R
│   │   ├── input.LOSO.rf.models.sh
│   │   ├── run_SHAP_test.R
│   │   └── SHAP_job_list.sh
│   ├── plotting/
│   │   ├── Extended.data.figure1.R
│   │   ├── Extended.data.figure2.R
│   │   ├── Extended.data.figure3.R
│   │   ├── Extended.data.figure4.R
│   │   ├── Extended.data.figure5.R
│   │   ├── Figure1.R
│   │   ├── Figure2.R
│   │   ├── Figure3.R
│   │   ├── Figure4.R
│   │   ├── Figure5.R
│   │   └── Figure6.R
│   ├── cancerness.utils.R
│   ├── parameters.yml
│   ├── requirements.R
│   └── utils.R
└── README.md
```

> **Note:** The data are currently only accessible from within the EMBL intranet, but will be made publicly available on **Zenodo** or another suitable repository soon.

## Workflow Overview

### 1. Model Training (Random Forest)

-   Machine learning models (Unified CRC, EO-CRC, LO-CRC, etc.) are trained using scripts in the `models/` directory.
-   Output files from model training are saved to `data/results/`.
-   To run LOSO models and evaluate both CRC and AD samples, use the scripts:
    -   `input.LOSO.rf.models.sh`
    -   `04.run.Train.LOSO.rf.models.sh`\
        These scripts must be submitted via SLURM.
-   To obtain study-to-study transfer (SST) evaluation AUCs, run:
    -   `05.Train.SCV.SST.rf.models.R`
-   Outputs from both LOSO and SST evaluations are stored in:\
    `data/results/scv.loso/`

### 2. Statistical Analysis & Linear Mixed Models

-   All statistical analysis scripts are located in `src/analysis/`.
-   Example scripts:
    -   `02.lmm.16S.WGS.R`: LMMs comparing CRC vs. CTR in 16S and WGS profiles
    -   `03.lmm.EO.LO.R`: LMMs for EO-CRC vs. EO-CTR or LO-CRC vs. LO-CTR comparisons

### 3. Model Evaluation & SHAP Analysis

-   SHAP interpretation scripts are located in `src/analysis/` and are submitted via SLURM (`04*.sh`).
-   Median SHAP values are computed and saved under:\
    `data/results/shap.analysis/`

### 4. Figure Generation

-   All manuscript figures are generated using scripts in `src/plotting/`, e.g.:
    -   `Figure1.R`, `Figure2.R`, ..., `Figure6.R`

## Contact

For questions or collaborations, please contact:

Selin Pekel — [s.pekel\@lumc.nl](mailto:s.pekel@lumc.nl)

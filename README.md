# Whole Brain Cellular Resolution Functional Network Properties of Seizure Susceptibility

## Project Overview
This project examines genotype-specific network changes in zebrafish with *scn1lab* mutations after PTZ exposure. Using calcium imaging and statistical modeling, we identify discriminatory metrics that differentiate seizure susceptibility.

**Reference Publication:** *Whole brain cellular resolution functional network properties of seizure susceptibility* (April 9, 2026)  
**DOI:** [https://doi.org/10.64898/2026.04.04.715761](https://doi.org/10.64898/2026.04.04.715761)

## Repository Contents & Data Availability
This repository contains the prepacked data along with the analysis and visualization code used to produce the 7 figures in the main paper ([10.26188/32061549](https://doi.org/10.26188/32061549)). A complete description and pseudocode of the core algorithmic steps can be found in the Methods and Supplementary sections of the main manuscript.

**Note:** Due to the large size of the original data and the extensive processing time required to prepare it, the raw data and original processing code are not hosted here. 
* The original calcium imaging data supporting the findings of this study are stored on Melbourne Mediaflux, and a sharing link is hosted on Melbourne Figshare (10.26188/32060955). The stored imaging data are password-protected; however, the password can be provided upon reasonable request by contacting the authors.
* Processed datasets are publicly available via Melbourne Figshare (10.26188/32060970). 
* All simulated data are accessible via Figshare (10.26188/32061360). 

The preprocessing pipeline is available here: [Scott Lab HPC Pipeline]([#](https://github.com/Scott-Lab-QBI/hpc_pipeline)) 

---

## 1. System Requirements

* **Operating Systems:** Windows 10/11, macOS Monterey 12.0, or Ubuntu 20.04 LTS. 
* **Tested Versions:** The software has been tested specifically on Windows 10.
* **Hardware:** A standard desktop computer is sufficient. No non-standard hardware is required. Processing times will vary based on CPU and RAM.
* **Software Dependencies:**
  * **Python:** Used for ROI Extraction & Registration (Suite2p for source localization/motion correction, and ANTs for 3D template warping via the Scott Lab HPC Pipeline).
  * **MATLAB (R2017a or newer):** Used for ΔF/F normalization, network generation, and statistical analysis. (Tested on version R2023b).
  * **R (version 4.4.1):** Used to produce specific plots and perform supplementary statistical analysis.
* **Required MATLAB Toolboxes:**
  * Statistics and Machine Learning Toolbox
  * Parallel Computing Toolbox (Required for generative modeling and simulation scripts)
  * Image Processing Toolbox
* **Required Third-Party Libraries:**
  * **Brain Connectivity Toolbox (BCT):** Crucial for graph theory and topological metrics. You must download it from the [BCT Website](https://sites.google.com/site/bctnet/).

## 2. Installation Guide

1. Clone this repository to your local machine:
   ```bash
   git clone [https://github.com/qinwayne/scn1lab.git](https://github.com/qinwayne/scn1lab.git)
   cd scn1lab

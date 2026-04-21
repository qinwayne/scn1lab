# scn1lab

### Project Overview
This project examines genotype specific network changes in zebrafish with scn1lab mutations after PTZ exposure. Using calcium imaging and statistical modeling, we identify discriminatory metrics that differentiate seizure susceptibility. 

**Reference Publication:**
Whole brain cellular resolution functional network properties of seizure susceptibility (April 9, 2026)
DOI: https://doi.org/10.64898/2026.04.04.715761

### Repository Contents
This repository contains the prepacked data along with the analysis and visualization code used to produce the 7 figures in the main paper.

Note: Due to the large size of the original data and the extensive processing time required to prepare it, the raw data and original processing code are not hosted here. The original calcium imaging data supporting the findings of this study are stored on Melbourne Mediaflux, and a sharing link is hosted on Melbourne Figshare 10.26188/32060955. Processed datasets are also publicly available via Melbourne Figshare 10.26188/32060970. All simulated data are accessible via Figshare 10.26188/32061360. The stored imaging data are password-protected; however, the password can be provided upon reasonable request by contacting the authors.

The preprocessing pipeline is available here: [Scott Lab HPC Pipeline](https://github.com/Scott%2DLab%2DQBI/hpc%2Dpipeline)

### 1. System Requirements

* **Operating Systems:** Windows 10, macOS Monterey 12.0, or Ubuntu 20.04 LTS.
* **Software Dependencies:** ROI Extraction & Registration: Custom Python scripts utilizing Suite2p (for source localization and motion correction) and Advanced Normalization Tools (ANTs) (for 3D template warping). Signal Normalization: MATLAB (version R2023b) was used for ${\Delta}F/F$ normalization and subsequent statistical analysis. R (version 4.41) was also used to produce plots and perform statistical analysis.
* **Tested Versions:** The software has been tested on Windows 11.
* **Hardware:** A standard desktop computer is sufficient but the required processing time varies. No non standard hardware is required.

### 2. Demo

* **Instructions to Run:** Execute the main visualization script provided in the folder which will automatically load the prepacked data.
* **Expected Output:** The code will output the 7 main figures presented in the publication, saved as high resolution image files.
* **Run Time:** The expected run time to generate all figures on a normal desktop computer is approximately 5 to 15 minutes.

* **License:** This software is provided under the MIT License, which is approved by the Open Source Initiative.
* **Repository Link:** https://github.com/qinwayne/scn1lab
* **Code Functionality:** A complete description and pseudocode of the core algorithmic steps can be found in the Methods section of the main manuscript.

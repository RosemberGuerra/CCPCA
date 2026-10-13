# CCPCA: Comparison of Penalized and Constrained Sparse PCA

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue.svg)](https://www.r-project.org/)
[![DOI](https://img.shields.io/badge/DOI-10.1007%2Fs11634--022--00499--2-blue)](https://doi.org/10.1007/s11634-022-00499-2)

> **Replication repository for simulation studies and empirical comparisons between penalized and constrained Sparse PCA formulations.**

This repository contains reproducible R scripts to replicate the numerical simulations and empirical applications published in:

> Guerra-Urzola, R., Van Deun, K., Vera, J. C., & Sijtsma, K. (2022). **A comparison of penalized and constrained approaches to sparse principal component analysis**. *Advances in Data Analysis and Classification*, 17(2), 411–442. [https://doi.org/10.1007/s11634-022-00499-2](https://doi.org/10.1007/s11634-022-00499-2)

---

## Repository Structure

```text
├── Application/
│   ├── Application.R          # Real-world data application (Section 4)
│   └── ccpca_functions.R      # Helper functions for the empirical application
├── Simulation/
│   ├── datageneration.R       # Generates simulation configurations across factor levels
│   ├── makeData.R             # Low-level synthetic data generator
│   ├── estimation.R           # Core benchmark comparing penalized vs. constrained SPCA
│   ├── analysing_results.R    # Summary evaluation and performance metrics
│   └── ccpca_functions.R      # Core algorithm routines and optimization solvers
└── README.md
```
## Getting Started & Prerequisites
Ensure you have R (>= 4.0.0) installed along with any required packages. Clone this repository locally:
```bash
git clone https://github.com/RosemberGuerra/CCPCA.git
cd CCPCA
```
(Note: The scripts automatically generate intermediate directories Data/ and Output/ if configured, or you can create them before running: mkdir -p Simulation/Data Simulation/Output.)

## Reproduction Workflow
### 1. Simulation Study (Section 3)
To replicate the simulation benchmarks comparing model formulations across noise levels and dimensions:

Generate Synthetic Datasets:

```R
# Generates simulated blocks across experimental conditions
source("Simulation/datageneration.R")
```
Datasets are written to the Simulation/Data/ directory.

Run Estimations:

```R
# Fits penalized and constrained Sparse PCA models
source("Simulation/estimation.R")
```
Model outputs are saved in Simulation/Output/.

Analyze & Summarize Results:

```R
# Computes evaluation metrics and generates comparative summary statistics
source("Simulation/analysing_results.R")
```
2. Empirical Application (Section 4)
To replicate the real-world application study:

```R
source("Application/Application.R")
```

### 2. Empirical Application (Section 4)
To replicate the real-world application study:

```R
source("Application/Application.R")
```
## Citation
If you use these algorithms or simulation benchmarks in your research, please cite the published paper:

```
@article{guerra2023sparsifying,
  title={Sparsifying the least-squares approach to PCA: comparison of lasso and cardinality constraint},
  author={Guerra-Urzola, Rosember and de Schipper, Niek C and Tonne, Anya and Sijtsma, Klaas and Vera, Juan C and Van Deun, Katrijn},
  journal={Advances in Data Analysis and Classification},
  volume={17},
  number={1},
  pages={269--286},
  year={2023},
  publisher={Springer}
}
```
## License
This repository is distributed under the MIT License - see the LICENSE file for details.
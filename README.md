# CCPCA
This repository contains the codes needed for results replication of: Comparison of Penalized and Constrained Sparse PCA (DOI
https://doi.org/10.1007/s11634-022-00499-2).

## Simulation
The scrips in simulation folder are used to replicate the simulation study (Sect. 3 in the manuscript).
### Data Generation
The file "datageneration.R" should be runnned first. This files creates data sets with different setting for the factor and leves. "datageneration.R" makes used of "makeData.R"
The data is saved in a Data folder.
### Estimation 
The script "estimation.R" is used for the estimation. The results are saved in a Output folder.
The analysis of the estimation is carried out using "analysing_results.R"
## Application
To replacate the application (Sect. 4 in the manuscript), the script "Application.R" can be runned.

# White Matter Volume and Microstructural Integrity Associated with Fatigue in Relapsing Multiple Sclerosis  

This repository contains the data and R scripts used in the analyses presented in the paper:  

**Alejandra Figueroa-Vargas A., Navarrete S., et al.**  
**White Matter Volume and Microstructural Integrity Are Associated with Fatigue in Relapsing Multiple Sclerosis.**  
Under Review in *Scientific Reports*. 
Preprint available at: [https://www.researchsquare.com/article/rs-5419035/v1](https://www.researchsquare.com/article/rs-5419035/v1)  
 

## Repository Structure  

The repository is organized into the following folders:  

### `MRI/`  
This folder includes the R scripts for analyses based on data preprocessed with FreeSurfer and MRtrix3.  

#### Scripts  
- **`HDIofMCMC.r`**: Calculates the High-Density Interval (HDI) for the posterior distributions of Markov chains.  
- **`lasso_FA-fatiga.R`**: Implements a LASSO model to explore the relationship between fatigue and fractional anisotropy (FA) in patients.  
- **`lasso_interc_cognitiva.R`**: Fits a LASSO model between MRI parameters and cognitive fatigue, controlling for physical fatigue.  
- **`lasso_interc_fisica.R`**: Fits a LASSO model between MRI parameters and physical fatigue.  
- **`lasso_interc.R`**: Fits a LASSO model between MRI parameters and total fatigue.  
- **`PCA_grosor_v3.R`**: Computes principal components of MRI data.  

### `DATA/`  
This folder contains the following data files:  
- Precalculated FreeSurfer data:  
  - `data_seg_GM_WM_test.csv`  
- Precalculated MRtrix3 data:  
  - `Pacientes_ROI_data.xlsx`  
- Test scores:  
  - Includes data for MFIS, GAD7, PHQ9, PASAT, SDMT, EDSS, and medication usage (`Datos_paper.csv`).  
- Results from preprocessing scripts:  
  - `PCA_sorted.RData`  
  - `pca_20_components_grupal_wm.RData`  
  - `pca_20_components_grupal_gm.RData`  
  - `pca_10_components_grupal_aseg.RData`  

## How to Use  

1. Clone this repository:  
   ```bash
   git clone https://github.com/neurocics/Figueroa-Vargas_Navarrete_2025_Scientific_Reports.git

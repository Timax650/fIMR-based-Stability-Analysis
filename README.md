# fIMR-based Stability Analysis
Here are models and codes used in our paper "Frequency-spectra Impedance Margin Ratio based Stability Analysis of IBR-Penetrated Systems" for future publication.
The codes are revised from Simplus Grid Tool: [Simplus Grid Tool](https://github.com/Future-Power-Networks/Simplus-Grid-Tool).

## Installation
Run InstallSimplusGT.m to install Simplus Grid Tool.

## Models
Data and Simulink models in our paper can be found in 'Simplus-Grid-Tool\IMRF_case'.
1) IEEE68_bad: the original case.
2) IEEE68_bad_revised_para: the case after revising parameters.

## Calculate whole sytem admittances and fIMR curves
1) Open Matlab, and include all folders and subfolders of this repo in Matlab's path.
2) Open UserMain.m in the root path of the repo, select the case, and click run to plot curves.

## Stability analysis
1) Open CriticalIMRF.m in the root path and select the mode frequencies and the apparatus of interest.
2) Click run to plot contribution factors, sensitivities and the heat map.
# Code & Data for multiplex wormhole publication

## Contents:
### Code:
* **Section3_Assess_existing_panels.py**: Python code for assessing dimer loads, thermodynamics of off-target amplification, in existing panels (Section 3 in the paper).
* **Section3B_Existing_panel_offtarget_amp.R**: R code for predicting, filtering, and plotting amplification of non-target species (diet items, cooccurring carnivores) and off-target amplification within the target species genome (Section 3 in the paper). 
* **Section4A_Size_Simulations.py**: Python code for running simulations for varying input dataset sizes, output panel sizes, and algorithms (Section 4 in the paper).
* **Section4B_MultiLevel_Simulations.py**: Python code for running multi-level optimization simulations (Section 4 in the paper).
* **MainPlots.R**: Code used to produce figures 2, 3, 4, 6, 7 for the paper.

### Data for Simulations:
* Marten_Templates.csv
* Fisher_Templates.csv
* GrayFox_Templates.csv
* SeaOtter_Templates.csv

### Simulation Results:
* **MW_SIMULATIONS_MERGED.csv**: Dimer load summaries from varying input dataset size/output panel size/algorithm simulations.
* **ORDER_SUMMARY_ALL.csv**: Dimer load summaries from multi-level optimizations.

### Panel Assessment Results:
* SNP_panel_assessments.csv: Predicted dimer loads for existing panels, and predicted dimer loads after removing top 20% of offending primer pairs (Section 3 in the paper).
* Predicted_Primer_Binding.zip: Folder containing full runPrimerTree outputs for existing panels.


## Citation
multiplex wormhole: a pipeline for in silico design and evaluation of genotyping-by-amplicon sequencing panels. (*in prep*)

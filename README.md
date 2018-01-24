
README for the sulcal pits analysis performed in
Le Guen et al,
Genetic Influence on the Sulcal Pits: On the Origin of the First Cortical Folds
Cerebral Cortex, https://doi.org/10.1093/cercor/bhx098
Published: 21 April 2017

Execute scripts in the following order:

pits_DPF_distribution.py
Gather in .csv the DPF distribution of all the pits in each areal

plot_DPF_distribution.py
Plot histograms presenting the DPF distribution in each areal
Select the best threshold per areal using a mixture of 1 or 2 gaussians

pits_extraction.py
Create the phenotype files using (or not) the previous thresholds

run_solar.py
First, create the SOLAR PEDIGREE
Second, add the covariates in the phenotype dataframe
Then, perform SOLAR heritability analysis
Finally, parse the SOLAR h2 and pval estimates in python dictionaries

map_pits_heritability.py
Display the heritability estimate per hemisphere on Freesurfer fsaverage_sym,
as well the associated p-values.

makeped.tcl
load pedigree with SOLAR

pheno_analysis.tcl
Set the appropriate covariates and performed heritability analysis with SOLAR

### Estimating time of HIV infection from NGS data
by Vadim Puller, Richard Neher, and Jan Albert, biorxiv, [doi:10.1101/129387](https://doi.org/10.1101/129387)

Code associated with our manuscript on estimating the time of HIV infection from next generation sequencing data. 
This project uses whole genome deep sequencing data from [Zanini et al 2015](http://hiv.biozentrum.unibas.ch) to establish a method to estimate time of infection from viral diversity.


`Frequency_data`
[//]: # (The directory `Frequency_data` contains arrays of nucleotide frequencies in each of the samples on which the analysis is based.)
The directory `Frequency_data` contains the patient data on which the analysis is based.
These file are derived from the NGS data by mapping. 

`patient_tt.npy` :   time points for patient samples

`patient_data.npy` and `patient_mask.npy` : data and mask for the array of nucleotide frequencies for all time points

`patient_viral_load.npy` :  viral load data

`patient_dilutions.npy` and `patient_dilutions_mask.npy` : dilutions data (see [Zanini et al 2015])

`annotations.txt` : auxhiliary gene annotations info


The scripts `EDI_functions.py`, `EDI_plotting.py`, and `EDI_median_regression.py` contain the source code to load the data, determine the regression coefficients, and generate the figures in the manuscript.
`EDI_functions.py` : functions generating data from 
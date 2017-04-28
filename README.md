### Estimating time of HIV infection from NGS data
by Vadim Puller, Richard Neher, and Jan Albert, biorxiv, [doi:10.1101/129387](https://doi.org/10.1101/129387)

Code associated with our manuscript on estimating the time of HIV infection from next generation sequencing data. 
This project uses whole genome deep sequencing data from [Zanini et al 2015](http://hiv.biozentrum.unibas.ch) to establish a method to estimate time of infection from viral diversity.

`Frequency_data`
The directory `Frequency_data` contains arrays of nucleotide frequencies in each of the samples on which the analysis is based. 
These file are derived from the NGS data by mapping. 

The scripts `EDI_functions.py`, `EDI_plotting.py`, and `EDI_median_regression.py` contain the source code to load the data, determine the regression coefficients, and generate the figures in the manuscript.
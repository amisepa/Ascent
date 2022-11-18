# get_entropy()

EEGLAB plugin to compute different types of entropy measures from electrophysiology time series (e.g., EEG, MEG, ECG, HRV).

## Entropy measures

- Approximate entropy
- Sample entropy
- Fuzzy entropy
- Multiscale entropy
- Multiscale fuzzy entropy
- Refined composite multiscale fuzzy entropy (default; Azami and Escudero, 2017)

## Requirements

MATLAB, EEGLAB, some biosignal data

## Features

All parameters (tau, embedding dimension, fuzzy power, EEG channels of interest) can be adjusted easily in the graphical user interface (GUI) or in command line. For multiscale entropies, the mean, standard deviation, and variance can be used during the coarse-grain process. 

The option to control for spectral bias using bandpass filters at each time scale is available (see Kosciessa et al. 2020 for more detail).

Users can visualize the entropy indexes for each file, displayed on a scalp topography plot. For multiscale entropies, users can click on the sensor of interest to visualize the entropy at each time scale for that sensor. 

When possible, the code implements Matlab's parrallel computing to enhance computing speed.

## Tutorial

1) Install EEGLAB and the *get_entropy()* plugin

2) Tools > Compute entroy
<img src="https://github.com/amisepa/get_entropy/tree/main/tutorial/img2.png" width="400">

Here, select entropy measure of interest, EEG channels of interest, *tau*, *m*, and if you want to plot the result. 

Note that you need to import channel locations first if you want to plot the outputs (Edit > Channel locations). 


## References

Costa, Goldberger, and Peng (2005). Multiscale entropy analysis of biological signals. Phys. Rev. E 71, 021906.

Azami and Escudero (2017). Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis. Medical & Biological Engineering & Computing.

Kosciessa, Kloosterman, & Garrett (2020). Standard multiscale entropy reflects neural dynamics at mismatched temporal scales: What’s signal irregularity got to do with it? PLoS computational biology, 16(5), e1007885.

## Copyright 

Cedric Cannard, 2022

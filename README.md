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

2) Import your data and prepross them if necessary (import channel locations if you want to plot the entropy outputs).

3) Tools > Compute entroy

<img src="https://github.com/amisepa/get_entropy/blob/main/tutorial/img2.png" width="300">

4) 1st GUI window

<img src="https://github.com/amisepa/get_entropy/blob/main/tutorial/img3.png" width="300">

Here, you can select the entropy measure you wish to compute, the EEG channels of interest, *tau*, *m*, and if you want to plot the result. 

Note that you need to import channel locations before this step if you want to plot the outputs (Edit > Channel locations). 

5) If you selected a fuzzy entropy measure, this window will pop-up to allow you to adjjust fuzzy power. 

<img src="https://github.com/amisepa/get_entropy/blob/main/tutorial/img4.png" width="150">

6) If you selected a measure with multiple scale factors, a new GUI window opens:

<img src="https://github.com/amisepa/get_entropy/blob/main/tutorial/img6.png" width="300">

Here you can select the coarse graining method, the number of time scales, and if you want to control for spectral bias (i.e. apply bandpass filter to each scale; see Kosciessa et al., 2020). 

7) Plot for entropy measures with one value per channel (i.e., approximate, sample, and fuzzy entropy):

<img src="https://github.com/amisepa/get_entropy/blob/main/tutorial/img5.png" width="300">

You can rotate the 3D image if some channels are hard to see.

https://user-images.githubusercontent.com/58382227/202803880-f266aa24-77b8-449f-89ab-66d74e29b092.mp4


8) Plot for entropy measures with several values per channel (i.e., multiple time scales):

You can rotate the 3D image if some channels are hard to see and click on the channels to display the entropy values at all time scales.

https://user-images.githubusercontent.com/58382227/202819448-dc521f09-1229-4271-8eb3-662214597fec.mp4



## References

Costa, Goldberger, and Peng (2005). Multiscale entropy analysis of biological signals. Phys. Rev. E 71, 021906.

Azami and Escudero (2017). Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis. Medical & Biological Engineering & Computing.

Kosciessa, Kloosterman, & Garrett (2020). Standard multiscale entropy reflects neural dynamics at mismatched temporal scales: What’s signal irregularity got to do with it? PLoS computational biology, 16(5), e1007885.

## Copyright 

Cedric Cannard, 2022

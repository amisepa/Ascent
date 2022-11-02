# get_entropy()

EEGLAB plugin to compute different types of entropy measures:
- Approximate entropy
- Sample entropy
- Fuzzy entropy
- Multiscale entropy
- Multiscale fuzzy entropy
- Refined composite multiscale fuzzy entropy (default; Azami and Escudero, 2017)

All parameters (tau, embedding dimension, fuzzy power, EEG channels of interest) can be adjusted easily in the graphical user interface (GUI) or in command line. For multiscale entropies, the mean, standard deviation, and variance can be used during the coarse-grain process. 

The option to control for spectral bias using bandpass filters at each time scale is available as well (see Kosciessa et al. 2020 for more detail).

Users can visualize the entropy indexes for each file, displayed on a scalp topography plot. For multiscale entropies, users can click on the sensor of interest to visualize the entropy at each time scale for that sensor. 

Copyright - Cedric Cannard, 2022

# ASCENT EEGLAB plugin: Aperiodic Signal Complexity Estimation for Neurophysiological Time series

The Escape EEGLAB plugin computes entropy and complexity measures from multidimensional M/EEG data (or any other time series).

## Entropy measures available

Uniscale measures:
- Sample entropy (SampEn)
- Extrema-segmented entropy analysis of time series (ExSEnt; improved SampEn from Kamali 2025)
- Fuzzy entropy (FuzzEn)
- Fractal Dimension/Volatility (FracDim)

Multiscale measures:
- Multiscale entropy (MSE; enhanced version of Costa 2002)
- Modified Multiscale entropy (mMSE; enhanced version of Kloosterman 2020 and Kosciessa 2020)
- Multiscale fuzzy entropy (MFE; enhanced version of Azami 2017)
- Refined composite multiscale fuzzy entropy (RCMFE; enhanced version of Azami 2017)


## Requirements

- MATLAB
- EEGLAB


## Please cite 

https://psyarxiv.com/xwmyk/

## Illustration 

Here, we computed all available complexity measures from two conditions of 64-channel Biosemi data: eyes-open vs eyes-closed resting state.
We then perform (1,000 iterations) bootstrap statistics to identify the significant differences under the null hypothesis (H0; α = 0.05), and apply threshold-free cluster enhancement (TFCE) correction to control for the family-wise-error (FWE; i.e. Type 1 error), highlighting the significant spatioemporal clusters. 
For the multiscale measures, This for two coarsing methods (standard deviation and median)  to clearly visualize how they capture different types of complexity information. 
Time to compute everything with 32 GB of RAM and 10 CPU cores (parallel computing ON): ~11 hours. 

### Uniscale measures

<img width="2310" height="1824" alt="uniscales" src="https://github.com/amisepa/Ascent/blob/main/figures/uniscales.png" />


### Multiscale measures

<img width="3245" height="1847" alt="multiscale_std" src="https://github.com/amisepa/Ascent/blob/main/figures/multiscale_std.png" />

<img width="3240" height="1727" alt="multiscale_median" src="https://github.com/amisepa/Ascent/blob/main/figures/multiscale_median.png" />


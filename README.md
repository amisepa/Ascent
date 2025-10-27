# ASCENT EEGLAB plugin: Aperiodic Signal Complexity Estimation for Neurophysiological Time series

The Escape EEGLAB plugin computes entropy and complexity measures from multidimensional M/EEG data (or any other time series).

## Entropy measures available

Uniscale measures:
- Sample entropy (SampEn)
- Fuzzy entropy (FuzzEn)
- Fractal Dimension/Volatility (FracDim)

Multiscale measures:
- Multiscale entropy (MSE; enhanced version of Costa 2002)
- Modified Multiscale entropy (mMSE; enhanced version of Kloosterman 2020 and Kosciessa 2020)
- Multiscale fuzzy entropy (MFE; enhanced version of Azami 2017)
- Extrema-segmented entropy analysis of time series (ExSEnt; enhanced version of Kamali 2025)
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

<img width="2310" height="1824" alt="uniscales" src="https://github.com/user-attachments/assets/fe66bcdc-a5c1-487d-94ae-ffc778de5567" />

### Multiscale measures

<img width="3245" height="2378" alt="multiscale_std" src="https://github.com/user-attachments/assets/00c40380-c46e-4bd6-b9db-7a0b8fae8f40" />

<img width="3245" height="2323" alt="multiscale_median" src="https://github.com/user-attachments/assets/2f3968c9-36f2-4adc-b2c0-b0673b22c3f6" />




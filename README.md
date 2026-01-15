
<p align="center" width="100%">
  <img width="50%" alt="BrainBeats logo"
       src="https://raw.githubusercontent.com/amisepa/Ascent/main/ascent_ai_img2.png">
</p>


The ASCENT EEGLAB plugin computes entropy and complexity measures from multidimensional M/EEG data (or any other time series).

## Measures supported

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

All algorithms were modified to significantly increase computation speed while preserving the integrity of the complexity estimates via vectorization, matrix operations, parallel computing, and blockwise bounded-memory blockwise distance calculations to avoid full pairwise matrix allocation. 


## Requirements

- MATLAB
- EEGLAB


## Please cite 

https://psyarxiv.com/xwmyk/

## Illustration 

We computed all of ASCENT's measures for two conditions of 64-channel Biosemi data: eyes-open vs eyes-closed resting state (N = 40).
We then performed (5,000 iterations) bootstrap statistics to identify the significant differences under the null hypothesis (H0; α = 0.05), and applied threshold-free cluster enhancement (TFCE) correction to control for the family-wise-error (FWE; Type 1 error), highlighting the significant spatiotemporal clusters. 

For the multiscale measures, group analysis was performed using the mean as the coarse-graining method. 

Time to compute everything with 32 GB of RAM and 10 cores with parallel computing: ~11 hours. 

### Uniscale measures

<img width="2310" height="1824" alt="uniscales" src="https://github.com/amisepa/Ascent/blob/main/figures/fig2.png" />

<img width="3245" height="1847" alt="multiscale_std" src="https://github.com/amisepa/Ascent/blob/main/figures/fig3.png" />


### Multiscale measures

<img width="3240" height="1727" alt="multiscale_median" src="https://github.com/amisepa/Ascent/blob/main/figures/fig4.png" />



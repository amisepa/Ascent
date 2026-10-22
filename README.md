
<p align="center" width="100%">
  <img width="30%" alt="ASCENT logo" src="img/ascent_logo.png">
</p>

# ASCENT

**Aperiodic Signal Complexity Estimation for Neurophysiological Time series**

ASCENT is an EEGLAB plugin to compute entropy, fractal, and aperiodic (1/f) measures on M/EEG data, either on the scalp channels or on the independent components (ICs). It can also be used on any other multichannel signal that EEGLAB can import (e.g., iEEG, ECoG, local field potentials). All measures can be computed from a graphical interface or from the command line, with the same parameters.

## Measures supported

Uniscale measures:
- Sample entropy (SampEn)
- Fuzzy entropy (FuzzEn)
- Extrema-segmented entropy (ExSEnt; improved SampEn from Kamali 2025)
- Fractal dimension/volatility (FracDim)
- Higuchi fractal dimension (HigFracDim)
- Aperiodic exponent and offset of the power spectrum (specparam algorithm from Donoghue 2020), static or time-resolved

Multiscale measures:
- Multiscale entropy (MSE; enhanced version of Costa 2002)
- Modified multiscale entropy (mMSE; enhanced version of Kloosterman 2020 and Kosciessa 2020)
- Multiscale fuzzy entropy (MFE; enhanced version of Azami 2017)
- Refined composite multiscale fuzzy entropy (RCMFE; enhanced version of Azami 2017)

Multivariate measures:
- Multivariate fuzzy entropy (mvFuzzEn; Ahmed 2011, Azami 2017)
- Refined composite multivariate multiscale fuzzy entropy (RCmvMFE; Azami 2017)

All algorithms were modified to significantly increase computation speed while preserving the integrity of the estimates, using vectorization, matrix operations, parallel computing, and blockwise distance calculations that avoid allocating the full pairwise matrix when memory is limited.

## Installation

Requires MATLAB and EEGLAB. Download or clone this repository into the `plugins` folder of EEGLAB and restart EEGLAB. An **ASCENT** menu then appears in the EEGLAB menu bar.

## Getting started

### Graphical interface

1. Start EEGLAB and click **ASCENT** in the menu bar.

<img alt="ASCENT in the EEGLAB menu bar" src="img/gui_tab.png" width="45%" />

2. If no dataset is loaded in EEGLAB, load one from the ASCENT window (*Load dataset...*), or click *Sample dataset* to load the 64-channel BioSemi recording provided with the plugin (resting state with eyes closed, with an ICA decomposition and dipole fits).

<img alt="ASCENT main window" src="img/gui_main.png" width="60%" />

3. Select the measure, the data to use (scalp channels or ICs), the channels (all by default), and whether to plot the outputs, track progress in the command window, and use parallel computing. Then click *Ok*.

4. A second window sets the parameters of the selected measure. Default values follow published guidelines (e.g., m = 2, r = 0.15). For the multiscale measures, users can also choose the coarse-graining method and the number of scale factors.

<img alt="Multiscale parameters" src="img/gui_multiscale.png" width="60%" />

For the aperiodic measure, this window sets the power spectrum estimation and the spectral parameterization (fixed or knee mode, number and width of the peaks), whether to remove the aperiodic component from the spectrum, and whether to also estimate it over time.

<img alt="Aperiodic parameters" src="img/gui_aperiodic.png" width="60%" />

5. The outputs are plotted and stored in `EEG.ascent.<measure>` (e.g., `EEG.ascent.RCMFE.data`, channels x scales). The dataset with its outputs is added to EEGLAB, and the command is saved in the EEGLAB history (type `eegh`), so that the same analysis can be run from a script. If the dataset was not modified since it was loaded, the outputs are also saved in its .set file.

For the multiscale measures, the figure shows the entropy of each channel at each scale with the high-entropy cluster outlined, the scalp topography of this cluster, and its mean entropy across scales. Below, RCMFE computed on the sample dataset:

<img alt="RCMFE on the sample dataset" src="img/example_rcmfe.png" width="70%" />

## Tutorial: each measure on the sample dataset

The commands below compute each measure on the sample dataset provided with the plugin (64-channel BioSemi, 75 s of resting state with eyes closed, with an ICA decomposition and dipole fits), with the settings used in the paper. Each call plots its outputs (`'vis', false` to turn this off) and stores them in `EEG.ascent.<measure>`. See `help ascent_compute` for all parameters and `ascent_tutorial.m` for more options.

```matlab
eeglab nogui
pluginPath = fileparts(which('eegplugin_ascent.m'));
EEG = pop_loadset('filename', 'ascent_sample_data.set', 'filepath', pluginPath);
```

### Single-scale measures

```matlab
EEG = ascent_compute(EEG, 'measure', 'SampEn');      % sample entropy
EEG = ascent_compute(EEG, 'measure', 'FuzzEn');      % fuzzy entropy
EEG = ascent_compute(EEG, 'measure', 'ExSEnt');      % extrema-segmented entropy (durations, amplitudes, joint)
EEG = ascent_compute(EEG, 'measure', 'FracDim');     % fractal dimension (box-counting)
EEG = ascent_compute(EEG, 'measure', 'HigFracDim');  % Higuchi fractal dimension

EEG.ascent.SampEn.data                               % one value per channel
```

Parameters can be tuned with name-value pairs, for example `'tau', 1, 'm', 2, 'r', 0.15` for the entropy measures, or `'chanlist', {'Fz' 'Cz' 'Pz'}` to select channels.

### Aperiodic activity (1/f)

```matlab
% Exponent and offset of the power spectrum (specparam algorithm), alpha power before
% and after removing the aperiodic component, and the same over time (sliding windows)
EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'timeResolved', true);

EEG.ascent.Aperiodic.data.exponent                   % one value per channel
EEG.ascent.Aperiodic.data.offset
EEG.ascent.Aperiodic.data.exp_slid                   % channels x time windows
```

Left: exponent and offset topographies, their correlation across channels, raw and aperiodic-corrected PSD, and alpha power before and after correction. Right: time-resolved spectrograms, exponent, and offset.

<img alt="Aperiodic on channels" src="img/fig4_aperiodic_channels.png" width="85%" />

### Multiscale measures

```matlab
EEG = ascent_compute(EEG, 'measure', 'MSE',   'coarsing', 'mean', 'num_scales', 30);  % multiscale entropy
EEG = ascent_compute(EEG, 'measure', 'mMSE',  'coarsing', 'mean', 'num_scales', 30);  % modified MSE
EEG = ascent_compute(EEG, 'measure', 'MFE',   'coarsing', 'mean', 'num_scales', 30);  % multiscale fuzzy entropy

% Refined composite multiscale fuzzy entropy, with mean and standard deviation
% coarse-graining (scales 2-50)
EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'coarsing', 'mean', 'num_scales', 50);
EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'coarsing', 'std',  'num_scales', 50);

EEG.ascent.RCMFE.data                                % channels x scales
EEG.ascent.RCMFE.scales
```

The figure shows the entropy of each channel at each scale, with the high-entropy cluster outlined in orange (channels whose entropy exceeds the 75th percentile over the cluster's scales), its scalp topography, and its mean entropy across scales. Top: mean coarse-graining. Bottom: standard deviation coarse-graining.

<img alt="RCMFE on channels" src="img/fig2_rcmfe_channels.png" width="70%" />

The cluster threshold can be changed when plotting again, e.g. `ascent_plot(EEG.ascent.RCMFE.data, EEG.chanlocs, 'RCMFE', EEG.ascent.RCMFE.scales, 'ClusterThresh', 0.9)`.

### Independent components and dipoles

Any measure can be computed on the independent components (ICs) instead of the scalp channels with `'domain', 'ica'`. IC results are back-projected to the scalp for the topographies, and when the dataset has dipole fits (DIPFIT), a second figure shows the equivalent dipoles of the ICs on the MNI template and the entropy by brain region.

```matlab
EEG = ascent_compute(EEG, 'measure', 'RCMFE', 'domain', 'ica', 'coarsing', 'std', 'num_scales', 30);
EEG = ascent_compute(EEG, 'measure', 'Aperiodic', 'domain', 'ica', 'timeResolved', true);
```

<img alt="RCMFE on ICs" src="img/fig3_rcmfe_ics_dipoles.png" width="70%" />

<img alt="Aperiodic on ICs" src="img/fig5_aperiodic_ics.png" width="85%" />

### Multivariate measures

These are computed jointly over the selected channels and can be very long with many channels, so we recommend using them on a few channels, regions of interest, or components:

```matlab
EEG = ascent_compute(EEG, 'measure', 'mvFuzzEn', 'chanlist', {'Fz' 'Cz' 'Pz' 'Oz'});
EEG = ascent_compute(EEG, 'measure', 'RCmvMFE',  'chanlist', {'Fz' 'Cz' 'Pz' 'Oz'}, 'coarsing', 'std', 'num_scales', 10);
```

## Please cite

https://psyarxiv.com/xwmyk/

## Main results (eyes closed vs eyes open, N = 40)

In the paper, we computed these measures on 64-channel BioSemi recordings of 40 participants at rest with eyes closed (EC) and eyes open (EO), and compared the two conditions with paired permutation t-tests (2,000 permutations) and cluster-based correction for the family-wise error (p < 0.05). The script is provided (`ascent_group_analysis.m`).

Single-scale measures (t-values; channels outside significant clusters are masked):

<img alt="Single-scale measures" src="img/fig6_uniscale_group.png" width="85%" />

Multiscale measures, with standard deviation coarse-graining:

<img alt="Multiscale, SD coarse-graining" src="img/fig8_multiscale_sd_group.png" width="70%" />

and with mean coarse-graining:

<img alt="Multiscale, mean coarse-graining" src="img/fig9_multiscale_mean_group.png" width="70%" />

Power spectral density, raw (first) and after removing the aperiodic component (second):

<img alt="PSD, raw" src="img/fig10_psd_raw_group.png" width="70%" />

<img alt="PSD, aperiodic-corrected" src="img/fig11_psd_corrected_group.png" width="70%" />

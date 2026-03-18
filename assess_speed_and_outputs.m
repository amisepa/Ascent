%% Calculate speed improvements made in ASCENT algorithms compared to
% original and validate outputs' integrity.


%% Sample Entropy

clear; close all; clc % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
addpath(genpath(pluginPath))
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = EEG.data;
m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = true; % parallel computing
trackProg   = true; % track progress
parpool


% Original from Azami (2017)
tic
out1 = [];
for iChan = 1:EEG.nbchan
    out1(iChan,:) = compute_SampEn_azami(data(iChan,:), m, r);
    fprintf('  ch %3d/%3d: %.3f\n', iChan, EEG.nbchan, out1(iChan,:));
end
toc

% Ascent
tic
out2 = [];
out2 = compute_SampEn(data, 'm', m, 'r', r, ...
    'parallel', paraComp, 'Progress', trackProg);
toc

% Plot
figure('color','w'); 
nexttile; hold on
plot(out1, 'xr', 'LineWidth',2); 
plot(out2, 'ob', 'LineWidth',2, 'MarkerSize', 8); 
legend("original", "Ascent")
xlabel("EEG channels"); ylabel("Entropy")
title("Sample Entropy")

% Check outputs
format long g
difference = out1 - out2

%% Fuzzy Entropy

clear % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
addpath(genpath(pluginPath))
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = EEG.data;
m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = true; % parallel computing
trackProg   = true; % track progress


% original from Azami (2017)
tic
out1 = [];
for iChan = 1:EEG.nbchan
    out1(iChan,:) = compute_FuzzEn_azami(data(iChan,:), m, 'Exponential', [r n], 0, tau); % (ts, m, mf, rn, local,tau)
    fprintf('  ch %3d/%3d: %.3f\n', iChan, EEG.nbchan, out1(iChan,:));
end
toc

% Ascent
tic
out2 = compute_FuzzEn(data, 'm', m, 'tau', tau, 'n', n, 'r', r, ...
    'Kernel', 'exponential', 'BlockSize', 256, ...
    'Parallel', paraComp, 'Progress', trackProg);
toc

% Plot
nexttile; hold on
plot(out1, 'xr', 'LineWidth',2); 
plot(out2, 'ob', 'LineWidth',2, 'MarkerSize', 8); 
legend("original", "Ascent")
xlabel("EEG channels"); ylabel("Entropy")
title("Fuzzy Entropy")

% Check outputs
difference = out1 - out2

%% Multivariate fuzzy entropy

clear % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
addpath(genpath(pluginPath))
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = EEG.data;
m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = true; % parallel computing
trackProg   = true; % track progress

% % Ascent
% [entropy, phi_m, phi_m1] = compute_mvFuzzEn(data, 'm', m, 'tau', tau, 'n', n, 'r', r, ...
%     'Kernel', kernel_meth, 'BlockSize',blocksize, ...
%     'Parallel', paraComp, 'Progress', trackProg);


%% Extrema-Segmented Entropy (ExSEnt)

clear % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
addpath(genpath(pluginPath))
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = EEG.data;
m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = true; % parallel computing
trackProg   = true; % track progress

% original
tic
out1 = [];
for iChan = 1:EEG.nbchan
    % out1(iChan,:) = compute_FuzzEn_azami(data(iChan,:), m, 'Exponential', [r n], 0, tau); % (ts, m, mf, rn, local,tau)
    [HD, HA, out1(iChan,:)] = compute_ExSEnt_ori(zscore(data(iChan,:)),0.001,m,r);
    fprintf('  ch %3d/%3d: %.3f\n', iChan, EEG.nbchan, out1(iChan,:));
end
toc


% Ascent
tic
[HD, HA, out2] = compute_ExSEnt(data, 'm', m, 'r', r, ...
    'lambda', 0.001, 'Plot', false, ...
    'Parallel', paraComp, 'Progress', trackProg);
toc

% Plot
nexttile; hold on
plot(out1, 'xr', 'LineWidth',2); 
plot(out2, 'ob', 'LineWidth',2, 'MarkerSize', 8); 
legend("original", "Ascent")
xlabel("EEG channels"); ylabel("Entropy")
title("ExSEnt")

% Check outputs
difference = out1 - out2

%% Fractal Dimension (FracDim)

[entropy, SD, info] = compute_FracDim(data, ...
    'RejectBursts', true, 'WinFrac', 0.02, 'ZThresh', 6, ...
    'RobustFit', 'theilsen', ... %  'theilsen' (default) or 'ols'
    'ScaleTrimIQR', true, ...
    'Parallel', paraComp, 'Progress', trackProg);

%% Classic Multiscale Entropy (MSE) from Costa but enhanced

[entropy, scales] = compute_MSE(data, 'm', m, 'tau', tau, ...
    'coarsing', coarsing, 'num_scales', num_scales, ...
    'Parallel', paraComp, 'Progress', trackProg);

%% Modified Multiscale Entropy (mMSE; from Kloosterman and Kosciessa, implemented in Fieldtrip)

[entropy, scales, info] = compute_mMSE(data, 'm', m, 'tau', tau, 'r', r, ...
    'coarsing', coarsing, 'num_scales', num_scales, ...
    'Parallel', paraComp, 'Progress', trackProg, ...
    'filter_mode', filter_mode, 'fs', fs, ...
    'TimeWin', TimeWin, 'TimeStep', TimeStep); % for time-resolved version

%% Multiscale Fuzzy Entropy (MFE)
[entropy, scales] = compute_MFE(data, 'm', m, ...
    'tau', tau, 'r', r, 'coarsing', coarsing, 'num_scales', num_scales, ...
    'Parallel', paraComp, 'Progress', trackProg);

%% Refined Composite Multiscale Fuzzy Entropy (RCMFE)

[entropy, scales] = compute_RCMFE(data, 'm', m, 'tau', tau, ...
    'r', r, 'n', n, 'coarsing', coarsing, 'num_scales', num_scales, ...
    'Parallel', paraComp, 'Progress', trackProg);

%% Refined Composite Multivariate Multiscale Fuzzy Entropy (RCMFE)

% original (optimized by Cedric on Dec 17, 2025 otherwise algo
% takes forever or crashes)
% tic
% entropy1 = compute_RCmvMFE_ori(data, m, r, n, tau, num_scales, coarsing, ...
%     trackProg);
% toc

% Ascent
[entropy, scales] = compute_RCmvMFE(data,'m', m, 'tau', tau, 'r', r, ...
    'coarsing', coarsing, 'num_scales', num_scales, ...
    'Parallel', paraComp, 'Progress', trackProg);

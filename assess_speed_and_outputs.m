%% Calculate speed improvements made in ASCENT algorithms compared to
% original and validate outputs' integrity.
%
% Both sides run serially (paraComp = false), so the times compare the
% algorithms; set paraComp = true to time ASCENT on a parallel pool. Run from the
% MATLAB desktop (EEGLAB starts in every section) on an otherwise idle machine,
% with this folder as the current folder so that this copy of the plugin is tested.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).


%% Sample Entropy

clear; close all; clc % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
fprintf('Benchmarking the plugin in %s\n', pluginPath);
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = EEG.data;
m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = false; % parallel computing
trackProg   = true; % track progress


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
if usejava('desktop')
    figure('color','w');
    nexttile; hold on
    plot(out1, 'xr', 'LineWidth',2);
    plot(out2, 'ob', 'LineWidth',2, 'MarkerSize', 8);
    legend("original", "Ascent")
    xlabel("EEG channels"); ylabel("Entropy")
    title("Sample Entropy")
end

% Check outputs
format long g
difference = out1 - out2
fprintf('SampEn: max |difference| = %.3g\n', max(abs(difference(:))));

%% Fuzzy Entropy

clear % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = EEG.data;
m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = false; % parallel computing
trackProg   = true; % track progress


% Local (ASCENT's default, removes each vector's own baseline) and global
modes = {1 'local'; 0 'global'};
for iMode = 1:size(modes,1)

    % original from Azami (2017)
    tic
    out1 = [];
    for iChan = 1:EEG.nbchan
        out1(iChan,:) = compute_FuzzEn_azami(double(data(iChan,:)), m, 'Exponential', [r n], modes{iMode,1}, tau); % (ts, m, mf, rn, local,tau)
        fprintf('  ch %3d/%3d: %.3f\n', iChan, EEG.nbchan, out1(iChan,:));
    end
    toc

    % Ascent
    tic
    out2 = compute_FuzzEn(data, 'm', m, 'tau', tau, 'n', n, 'r', r, ...
        'Mode', modes{iMode,2}, 'Kernel', 'exponential', 'BlockSize', 256, ...
        'Parallel', paraComp, 'Progress', trackProg);
    toc

    % Plot
    if usejava('desktop')
        nexttile; hold on
        plot(out1, 'xr', 'LineWidth',2);
        plot(out2, 'ob', 'LineWidth',2, 'MarkerSize', 8);
        legend("original", "Ascent")
        xlabel("EEG channels"); ylabel("Entropy")
        title("Fuzzy Entropy (" + modes{iMode,2} + ")")
    end

    % Check outputs
    difference = out1 - out2
    fprintf('FuzzEn %s: max |difference| = %.3g\n', modes{iMode,2}, max(abs(difference(:))));
end

%% Multivariate fuzzy entropy

clear % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = EEG.data;
m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = false; % parallel computing
trackProg   = true; % track progress

% % Ascent
% [entropy, phi_m, phi_m1] = compute_mvFuzzEn(data, 'm', m, 'tau', tau, 'n', n, 'r', r, ...
%     'Kernel', kernel_meth, 'BlockSize',blocksize, ...
%     'Parallel', paraComp, 'Progress', trackProg);


% EntropyHub's MvFuzzEn pools all embedded vectors and allocates [Nt x Nt]
% matrices, so it only runs on few channels and short data
nCh = 4; nSamp = 2000;
entropy = MvFuzzEn(double(data(1:nCh, 1:nSamp)).', 'm', m*ones(nCh,1), ...
    'tau', tau*ones(nCh,1), 'r', [r n], 'Norm', true);



%% Extrema-Segmented Entropy (ExSEnt)

clear % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = EEG.data;
m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = false; % parallel computing
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
if usejava('desktop')
    nexttile; hold on
    plot(out1, 'xr', 'LineWidth',2);
    plot(out2, 'ob', 'LineWidth',2, 'MarkerSize', 8);
    legend("original", "Ascent")
    xlabel("EEG channels"); ylabel("Entropy")
    title("ExSEnt")
end

% Check outputs
difference = out1 - out2
fprintf('ExSEnt: max |difference| = %.3g\n', max(abs(difference(:))));

%% Fractal Dimension (FracDim)

% [entropy, SD, info] = compute_FracDim(data, ...
%     'RejectBursts', true, 'WinFrac', 0.02, 'ZThresh', 6, ...
%     'RobustFit', 'theilsen', ... %  'theilsen' (default) or 'ols'
%     'ScaleTrimIQR', true, ...
%     'Parallel', paraComp, 'Progress', trackProg);

%% Multiscale Entropy (MSE)
clear % clear every time to avoid potential memory bias
eeglab; close;
pluginPath = fileparts(which('eegplugin_ascent.m'));
cd(pluginPath)
EEG = pop_loadset('filename','ascent_sample_data.set','filepath',fullfile(pluginPath));
data = double(EEG.data);

m           = 2;    % embedding dimension
r           = .15;  % similarity bound
n           = 2;    % fuzzy power
tau         = 1;    % time lag
paraComp    = false; % parallel computing
trackProg   = true;  % track progress
num_scales  = 5;
coarsing    = 'mean';

% Reference: MSEn from EntropyHub (Flood & Grimm, 2021; compute_MSE_ori.m), with
% the tolerance rescaled by the SD of each coarse-grained series (RadNew = 1)
Mobj.Func = @SampEn;
Mobj.m    = m;
Mobj.r    = r;
tic
out1 = nan(EEG.nbchan, num_scales);
for iChan = 1:EEG.nbchan
    out1(iChan,:) = compute_MSE_ori(data(iChan,:), Mobj, ...
        'Scales', num_scales, 'Methodx', 'coarse', 'RadNew', 1);
    fprintf('  ch %g/%g\n', iChan, EEG.nbchan);
end
toc

% Ascent, same tolerance rule
tic
[out2, scales] = compute_MSE(data, 'm', m, 'r', r, 'tau', tau, ...
    'coarsing', coarsing, 'num_scales', num_scales, 'zNorm', 1, 'IncludeScale1', true, ...
    'Parallel', paraComp, 'Progress', trackProg);
toc

% Compare differences (per scale, over channels)
difference = out1 - out2;
fprintf('MSE vs EntropyHub, scales %s: max |difference| = %s\n', mat2str(scales), mat2str(max(abs(difference), [], 1), 3));
if usejava('desktop')
    figure('color','w','ToolBar','none');
    imagesc(difference); colorbar; colormap(parula);
    xlabel('Scale'); ylabel('Channel');
    title('MSE difference: EntropyHub minus ASCENT');
end


%% Modified Multiscale Entropy (mMSE; from Kloosterman and Kosciessa, implemented in Fieldtrip)

% [entropy, scales, info] = compute_mMSE(data, 'm', m, 'tau', tau, 'r', r, ...
%     'coarsing', coarsing, 'num_scales', num_scales, ...
%     'Parallel', paraComp, 'Progress', trackProg, ...
%     'filter_mode', filter_mode, 'fs', fs, ...
%     'TimeWin', TimeWin, 'TimeStep', TimeStep); % for time-resolved version

%% Multiscale Fuzzy Entropy (MFE)
[entropy, scales] = compute_MFE(data, 'm', m, ...
    'tau', tau, 'r', r, 'coarsing', coarsing, 'num_scales', num_scales, ...
    'Parallel', paraComp, 'Progress', trackProg);

%% Refined Composite Multiscale Fuzzy Entropy (RCMFE)

[entropy, scales] = compute_RCMFE(data, 'm', m, 'tau', tau, ...
    'r', r, 'n', n, 'coarsing', coarsing, 'num_scales', num_scales, ...
    'Parallel', paraComp, 'Progress', trackProg);

%% Refined Composite Multivariate Multiscale Fuzzy Entropy (RCmvMFE)

% original (optimized by Cedric on Dec 17, 2025 otherwise algo
% takes forever or crashes)
% tic
% entropy1 = compute_RCmvMFE_ori(data, m, r, n, tau, num_scales, coarsing, ...
%     trackProg);
% toc

% Ascent
% Five channels: the pooled m+1 embedding grows with channels x samples
tic
[entropy, scales] = compute_RCmvMFE(data(1:5,:),'m', m, 'tau', tau, 'r', r, ...
    'coarsing', coarsing, 'num_scales', num_scales, ...
    'Parallel', paraComp, 'Progress', trackProg);

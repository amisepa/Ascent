%% Computes multiscale entropy (MSE)
%
% INPUTS:
%   signal: univariate signal - a vector of size 1 x N (the number of sample points)
%   m: embedding dimension
%   r: threshold (it is usually equal to 0.15 of the standard deviation of a signal - 
%       because we normalize signals to have a standard deviation of 1, here, r is usually equal to 0.15)
%   tau: time lag (it is usually equal to 1)
%   coarseType: 'Mean', 'Standard deviation' (default), 'Variance'
%   nScales: the number of scale factors (default = 15)
%   filtData: bandpass filter each scale factor to control for spectral 
%       bias (1) or not (0; default)
%   fs: sample rate (Hz)
%
% OUTPUTS:
%   mse: entropy values for each scale factor
%
% Ref:
%   [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on
%   Standard Deviation for Biomedical Signal Analysis", Medical & Biological
%   Engineering & Computing, 2016.
%   [2] M. Costa, A. Goldberger, and C.K. Peng, "Multiscale Entropy Analysis
%   of Complex Physiologic Time Series", Physical review letters, vol. 89, no. 6, p. 068102, 2002.
%
% Cedric Cannard, 2022

function [mse, freqs] = compute_mse(signal, m, r, tau, coarseType, nScales, filtData, fs)

% Max scale factor cannot be greater than Nyquist frequency
nf = fs/2;
if nScales >= nf
    warning(["Scale factor cannot be as high as signal's Nyquist frequency. Lowering it to " num2str(nf-1) ]);
end

% Signal is centered and normalized to have SD = 1
signal = signal-mean(signal);
signal = signal./std(signal);

% Simplify SD coarse-graining name
if contains(lower(coarseType), 'standard')
    coarseType = 'SD';
end

% if vis
%     figure('color','w')
% end

mse = nan(1,nScales);
freqs = nan(2,nScales);
for iScale = 1:nScales
    
%     fprintf('Scale %d \n', iScale)
    % make copy of signal in case it is bandpass-filtered at each scale
    sig = signal;
    
    % scale factor bounds
    highBound = (1/iScale).*nf + .05*((1./iScale).*nf);
    lowBound = (1/(iScale+1)).*nf - .05*((1./(iScale+1)).*nf);
    % Display the corresponding frequencies for this scale
    freqs(:,iScale) = [round(lowBound,3) round(highBound,3) ];
    disp(['scale ' num2str(iScale) ': ' num2str(round(lowBound,3)) ' - ' num2str(round(highBound,3)) ' Hz']);

    % Bandpass filter outside these bounds to control for spectral bias (if selected)
    if filtData
        disp('Applying bandpasss-filter remove spectral bias from outside these bounds.')

        if iScale > 1
            if highBound-lowBound > .05*nf    % for broad passbands: Chebyshev Type I zero-phase bandpass
                [b,a] = cheby1(4,1,highBound/nf,'low');
                % figure; freqz(b,a);
                sig = filtfilt(b,a,sig);
                [b,a] = cheby1(4,1,lowBound/nf,'high');
                sig = filtfilt(b,a,sig);
            else                                 % Butterworth zero-phase bandpass
                [b,a] = butter(10,highBound/nf,'low');
                sig = filtfilt(b,a,sig);
                [b,a] = butter(10,lowBound/nf,'high');
                sig = filtfilt(b,a,sig);
            end
        else
            [b,a]= butter(10,lowBound/nf,'high');   % only Butterworth highpass for scale 1
            %         highcutoff = nf;
            sig = filtfilt(b,a,sig);
        end

%         % Visualize filter effect on the power spectrum
%         [psd,f] = pwelch(sig,[],[],[],fs);
%         plot(f,psd); hold on;
%         title(num2str(iScale)); legend([num2str(lowcutoff) '-' num2str(highcutoff)]);

    end


    switch coarseType

        case 'Mean'
%             disp('Selected coarse-graining method: Mean')
            y = reshape(sig(1:floor(length(sig)/iScale)*iScale),iScale,[]);
            x = mean(y,'omitnan');
            mse(:,iScale) = compute_se(x, m, r, tau);

        case 'SD'

%             disp('Selected coarse-graining method: Standard deviation')
            y = reshape(sig(1:floor(length(sig)/iScale)*iScale),iScale,[]);
            x = std(y,'omitnan');
            mse(:,iScale) = compute_se(x, m, r, tau);

        case 'Variance'

%             disp('Selected coarse-graining method: Variance')
            y = reshape(sig(1:floor(length(sig)/iScale)*iScale),iScale,[]);
            x = var(y,'omitnan');
            mse(:,iScale) = compute_se(x, m, r, tau);

    end

%     if vis
%         plot(1:nScales, mse); hold on;
%     end
end

% if vis
%     nans = isnan(mse(1,:));
%     if nans(1)
%         xticks(2:nScales); xticklabels(join(string(freqs(:,2:end)),1)); xtickangle(45)
%         xlim([2 nScales]);
%     else
%         xticks(1:nScales); xticklabels(join(string(freqs),1)); xtickangle(45)
%     end
% end

% Cite references
disp('Please cite: ')
switch coarseType
    case 'Mean'
        disp('Costa, Goldberger, and Peng (2002) - Multiscale entropy analysis of complex physiologic time series. Physical review letters.')
    case 'SD'
        disp('   [1] Costa, Goldberger, and Peng (2002) - Multiscale entropy analysis of complex physiologic time series. Physical review letters.')
        disp('   [2] Azami and Escudero (2016) - Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis. Medical & Biological Engineering & Computing')
    case 'Variance'
        disp('   [1] Costa, Goldberger, and Peng (2002) - Multiscale entropy analysis of complex physiologic time series. Physical review letters.')
        disp('   [2] Azami and Escudero (2016) - Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis. Medical & Biological Engineering & Computing')
end
if filtData
    disp('Bandpass filters were applied to each scale factor to control for spectral bias, following recommendations by: ');
    disp('   [3] Kosciessa, Kloosterman, and Garrett (2020) - Standard multiscale entropy reflects neural dynamics at mismatched temporal scales: What''s signal irregularity got to do with it? Plos Computational Biology.')
end

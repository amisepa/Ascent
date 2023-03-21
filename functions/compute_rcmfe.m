% Calculates the refined composite multiscale fuzzy entropy (RCMFE) 
% whose coarse-graining uses standard deviation.
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
%   n: fuzzy power
%   usegpu: use GPU computing (1), or not (0)
%
% Outputs:
%   rcmfe: entropy values for each scale factor
%   scales: lower and upper frequency bounds of each time scale
%
% 
% Ref:
%   [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on
%   Standard Deviation for Biomedical Signal Analysis", Medical & Biological
%   Engineering & Computing, 2016.
%
% Cedric Cannard, 2022

function [rcmfe, freqs] = compute_rcmfe(signal,m,r,tau,coarseType,nScales,filtData,fs,n,usegpu)

% if not srate provided somehow, try to interpolate
% if ~exist('fs', 'var')
% 	fs = (1/(EEG.times(end)) - EEG.times(1));
% 	warning(fprintf('No sample rate inputted, sample rate estimated: %g', fs))
%     errordlg('You need to input the sample rate.'); return
% end

% Simplify SD coarse-graining name
if contains(lower(coarseType), 'standard')
    coarseType = 'SD';
end

% max scale factor
nf = fs/2; %Nyquist frequency
if nScales >= nf
	warning("Scale factor cannot be as high as signal's Nyquist frequency. Lowering it to %g", nf-1);
end

% Signal is centered and normalized to standard deviation 1
disp('Normalizing signal to SD = 1.')
signal = signal - mean(signal);
signal = signal ./ std(signal);

rcmfe = nan(1,nScales);
freqs = nan(2,nScales);
for iScale = 2:nScales
    
    fprintf('Scale %g \n', iScale);

    temp_A = [];
    temp_B = [];
    
    % Apply filters following recommendations from Kosciessa et al. (2020)
    if filtData
        % LOWPASS only approach
    %     highcutoff = (1/iScale)*nf;
    %     if highcutoff == nf
    %         highcutoff = highcutoff-1;
    %     end
    %     [b,a] = butter(6,highcutoff/nf);
    %     signal = filtfilt(b,a,signal);
    %     disp(['scale ' num2str(iScale) ': ' num2str(0) ' - ' num2str(round(highcutoff,1)) ' Hz']);
                
        % BANDPASS approach
        highcutoff = (1/iScale).*nf + .05*((1./iScale).*nf);
        lowcutoff = (1/(iScale+1)).*nf - .05*((1./(iScale+1)).*nf);
        freqs(:,iScale) = [round(lowcutoff,3) round(highcutoff,3)];
        disp(['scale ' num2str(iScale) ': ' num2str(round(lowcutoff,3)) ' - ' num2str(round(highcutoff,3)) ' Hz']);
    %     freqs(iScale,:) = sprintf('scale %d: %2.5g - %2.5g', iScale, round(lowcutoff,1), round(highcutoff,1));
    %     disp(freqs(iScale,:));
    %     freqs(iScale,:) = ['Scale ' num2str(iScale) ': ' num2str(round(lowcutoff,1)) ' - ' num2str(round(highcutoff,1))];
        
        if iScale > 1
            if highcutoff-lowcutoff > .05*nf    % for broad passbands: Chebyshev Type I zero-phase bandpass 
                [b,a] = cheby1(4,1,highcutoff/nf,'low');
                % figure; freqz(b,a);
                signal = filtfilt(b,a,signal);
                [b,a] = cheby1(4,1,lowcutoff/nf,'high');
                signal = filtfilt(b,a,signal);
            else                                 % Butterworth zero-phase bandpass
                [b,a] = butter(10,highcutoff/nf,'low');
                signal = filtfilt(b,a,signal);
                [b,a] = butter(10,lowcutoff/nf,'high');
                signal = filtfilt(b,a,signal);
            end
        else
            [b,a] = butter(10,lowcutoff/nf,'high');   % only Butterworth highpass for scale 1
    %         highcutoff = nf;
            signal = filtfilt(b,a,signal);
        end
    
        % look at PSD
    %     [Pxx,f] = pwelch(signal,[],[],[],fs);
    %     figure; plot(f,Pxx); title(num2str(iScale)); legend([num2str(lowcutoff) '-' num2str(highcutoff)]);
    end

    for ii = 1:iScale

        sig = signal(ii:end);
    
        % Coarse-graining process 
        y = reshape(sig(1:floor(length(sig)/iScale)*iScale),iScale,[]);
        switch coarseType
            case 'Mean'
                sigCoarse = mean(y,'omitnan');
            case 'SD'
                sigCoarse = std(y,'omitnan');
            case 'Variance'
                sigCoarse = var(y,'omitnan');
        end
        
        % Because of filtering, the scale-wise SD decreases relative to the global scale-
        % invariant similarity bound r [29]. Thus, r is recalculated for each scale, 
        % thereby normalizing MSE with respect to changes in overall time series variation at each scale
        if filtData
            % r_adj = 0.5*std(sig);
            r_adj = r*std(sig);
        else
            r_adj = r;
        end

        % Compute fuzzy entropy on the coearsed signal
        [~, p] = compute_fe(sigCoarse,m,r_adj,n,tau);

        temp_A = [temp_A p(1)];
        temp_B = [temp_B p(2)];
        
    end
    
    % output
    A = sum(temp_A);
    B = sum(temp_B);
    rcmfe(iScale) = log(A/B);
    
end



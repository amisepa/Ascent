% Calculates the refined composite multiscale fuzzy entropy (RCMFE) 
% whose coarse-graining uses standard deviation (RCMFE_std)
% 
% Inputs:
%
% signal: univariate signal
% dim: embedding dimension
% thresh: threshold (default = 0.15 of signal's SD after normalization of signal to SD =1)
% fpower: fuzzy power (default = 2)
% tau: time lag (default = 1)
% nscales: number of scale factors (default = 30)
% srate = sample rate of signal;
%
% Outputs:
%
% Out_RCMFE: a vector showing the RCMFE_std of x
% freqs: frequencies corresponding to each time scale
%
% Ref:
% [1] H. Azami and J. Escudero, "Refined Multiscale Fuzzy Entropy based on Standard Deviation for Biomedical Signal Analysis", Medical & Biological Engineering &
% Computing, 2016.
%
% If you use the code, please make sure that you cite reference [1].
%
% Hamed Azami and Javier Escudero Rodriguez
% hamed.azami@ed.ac.uk and javier.escudero@ed.ac.uk
%
% Adapted by Cedric Cannard, 2022
% Changes to follow requirements by Kosciessa et al. (2020):
%   1) filtering out spectral contamination for each scale
%   2) adjust r to each scale

function [rcmfe, freqs] = get_rcmfe(signal,m,r,n,tau,nscales,fs)

% ADD INPUT CHECKS

% if not srate provided somehow, try to interpolate
if ~exist(fs, 'var')
	fs = (1/(EEG.times(end)) - EEG.times(1));
	warning([ 'No sample rate inputted, sample rate detected = ' num2str(fs) ])
end

% max scale factor
nf = fs/2; %Nyquist frequency
if nscales >= nf
	warning(["Scale factor cannot be as high as signal's Nyquist frequency. Lowering it to " num2str(nf-1) ]);
end

% Signal is centered and normalized to standard deviation 1
disp('Normalizing signal to SD = 1.')
signal = signal - mean(signal);
signal = signal ./ std(signal);

rcmfe = nan(1,nscales);
freqs = nan(2,nscales);

for iScale = 2:nscales
        
    temp_A = [];
    temp_B = [];
    
    % Apply filters following recommendations from Kosciessa et al. (2020)
    
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
        [b,a]= butter(10,lowcutoff/nf,'high');   % only Butterworth highpass for scale 1
%         highcutoff = nf;
        signal = filtfilt(b,a,signal);
    end

    %look at PSD
%     [Pxx,f] = pwelch(signal,[],[],[],fs);
%     figure; plot(f,Pxx); title(num2str(iScale)); legend([num2str(lowcutoff) '-' num2str(highcutoff)]);

    for ii = 1:iScale
        
        sig = signal(ii:end);
        
        % Coarse-graining process based on standard deviation
        y = reshape(sig(1:floor(length(sig)/iScale)*iScale),iScale,[]);
        sigCoarse = std(y);
        
        % Because of filtering, the scale-wise SD decreases relative to the global scale-
        % invariant similarity bound r [29]. Thus, r is recalculated for each scale, 
        % thereby normalizing MSE with respect to changes in overall time series variation at each scale
%         r_new = r*std(sig);
%         r_new = 0.5*std(sig);
        r_new = r*std(sig);

        % Compute fuzzy entropy
        [~,P] = FuzEn(sigCoarse,m,r_new,n,tau);
        temp_A = [temp_A P(1)];
        temp_B = [temp_B P(2)];
        
    end
    
    %output
    A = sum(temp_A);
    B = sum(temp_B);
    rcmfe(iScale) = log(A/B);
    
end



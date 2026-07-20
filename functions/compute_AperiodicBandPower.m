function [bp_raw, bp_osc, ap_model] = compute_AperiodicBandPower(freqs, psd, exponent, offset, band)
% Band power split into its total and oscillatory (aperiodic-removed) parts.
%
%   [bp_raw, bp_osc, ap_model] = compute_AperiodicBandPower(freqs, psd, ...
%       exponent, offset, [8 13])
%
% Single source of truth for the aperiodic ("1/f") model in Ascent. Anything
% that needs the model, a corrected spectrum, or a band power should call this
% rather than re-deriving 10.^(offset - exponent*log10(f)) locally — that
% formula was previously duplicated across ascent_compute, the sliding-window
% fit and the plotting layer, and the three had drifted apart.
%
% Inputs
%   freqs    : frequency vector [1 x nFreqs], Hz
%   psd      : power spectra [nSig x nFreqs], linear µV²/Hz (not log).
%              Rows are channels or ICs; matches EEGLAB convention.
%   exponent : aperiodic exponent per signal [nSig x 1], from compute_AperiodicFit
%   offset   : aperiodic offset per signal   [nSig x 1], log10 units
%   band     : [fmin fmax] Hz to average over (default [8 13], alpha)
%
% Outputs
%   bp_raw   : total band power per signal [nSig x 1], µV²/Hz.
%              Includes the aperiodic background, so it is dominated by
%              broadband power for signals with a steep 1/f — it is NOT a
%              measure of oscillatory activity on its own.
%   bp_osc   : oscillatory band power per signal [nSig x 1], µV²/Hz.
%              bp_raw minus the aperiodic model over the same band, i.e. the
%              power genuinely above the 1/f fit. Same unit as bp_raw, so the
%              two are directly comparable. Negative where the fit overshoots
%              the band (returned as-is; callers decide how to handle it).
%   ap_model : aperiodic model [nSig x nFreqs], linear µV²/Hz.
%              A corrected spectrum in the ratio sense is psd ./ ap_model.
%
% NOTE
%   Only valid for AperiodicMode = 'fixed'. The knee model is not a straight
%   line in log-log space and is not reconstructed here.
%
% Copyright (C) - ASCENT EEGLAB PLUGIN - Cedric Cannard, 2021-2025

if nargin < 5 || isempty(band), band = [8 13]; end

freqs    = freqs(:)';
exponent = exponent(:);
offset   = offset(:);

if size(psd,2) ~= numel(freqs)
    error('compute_AperiodicBandPower: psd has %d frequency bins but freqs has %d.', ...
          size(psd,2), numel(freqs));
end
if numel(exponent) ~= size(psd,1) || numel(offset) ~= size(psd,1)
    error('compute_AperiodicBandPower: exponent/offset must have one entry per psd row (%d).', ...
          size(psd,1));
end

% Aperiodic model in linear units: log10(P_ap) = offset - exponent*log10(f)
ap_model = 10.^(offset - exponent .* log10(freqs));   % implicit expansion -> [nSig x nFreqs]

mask = freqs >= band(1) & freqs <= band(2);
if ~any(mask)
    warning('compute_AperiodicBandPower: band [%g %g] Hz not covered by freqs [%g %g].', ...
            band(1), band(2), freqs(1), freqs(end));
    bp_raw = nan(size(psd,1),1);
    bp_osc = bp_raw;
    return
end

bp_raw = mean(psd(:,mask),      2, 'omitnan');
bp_osc = bp_raw - mean(ap_model(:,mask), 2, 'omitnan');
end

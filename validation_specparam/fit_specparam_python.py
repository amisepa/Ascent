# fit_specparam_python.py — specparam 2.0 reference fits for ASCENT validation
# Reads validation_specparam/matlab_results.mat (v7), fits the same spectra
# with matched settings, writes python_results.csv + python_full.mat.
import os, sys
import numpy as np
import scipy.io as sio
from specparam import SpectralModel

root   = os.path.expandvars(r'C:\Users\ccann\Documents\MATLAB\Ascent')
outdir = os.path.join(root, 'validation_specparam')
mat    = sio.loadmat(os.path.join(outdir, 'matlab_results.mat'), squeeze_me=True)

f_syn  = np.asarray(mat['f_ax'], dtype=float).ravel()          # 1:0.25:40
psdSyn = np.asarray(mat['psdSyn'], dtype=float)                # [nSyn x nF]
freqsR = np.asarray(mat['freqsR'], dtype=float).ravel()
psdReal = np.asarray(mat['psdReal'], dtype=float)              # [nChan x nF]

SETTINGS = dict(peak_width_limits=[1.0, 12.0], max_n_peaks=3,
                min_peak_height=0.05, peak_threshold=2.0)

def fit_batch(freqs, psds, tag):
    n = psds.shape[0]
    exp_ = np.full(n, np.nan); off_ = np.full(n, np.nan)
    r2_  = np.full(n, np.nan); mae_ = np.full(n, np.nan)
    npk_ = np.full(n, np.nan)
    for k in range(n):
        psd = psds[k]
        if not np.all(np.isfinite(psd)) or np.any(psd <= 0):
            continue
        try:
            m = SpectralModel(aperiodic_mode='fixed', verbose=False,
                              algorithm_settings=SETTINGS)
            m.fit(freqs, psd)
            p = m.results.params
            ap = np.asarray(p.aperiodic._fit, dtype=float).ravel()  # [offset, exponent]
            off_[k], exp_[k] = ap[0], ap[1]
            pe = np.asarray(p.periodic._fit, dtype=float)           # native gaussian [cf pw sigma]
            npk_[k] = pe.shape[0] if pe.size else 0
            mr = m.results.metrics.results          # {'error_mae': .., 'gof_rsquared': ..}
            mae_[k] = float(mr['error_mae']); r2_[k] = float(mr['gof_rsquared'])
        except Exception as e:
            print(f'{tag} {k}: FAIL {e}', flush=True)
        if (k+1) % 50 == 0:
            print(f'{tag}: {k+1}/{n}', flush=True)
    return exp_, off_, r2_, mae_, npk_

print('fitting synthetic...', flush=True)
eS, oS, r2S, maeS, npkS = fit_batch(f_syn, psdSyn, 'syn')
print('fitting real...', flush=True)
eR, oR, r2R, maeR, npkR = fit_batch(freqsR, psdReal, 'real')

sio.savemat(os.path.join(outdir, 'python_results.mat'), {
    'expPS': eS, 'offPS': oS, 'r2PS': r2S, 'maePS': maeS, 'npkPS': npkS,
    'expPR': eR, 'offPR': oR, 'r2PR': r2R, 'maePR': maeR, 'npkPR': npkR})
print('saved python_results.mat', flush=True)
print('==== PYTHON SIDE DONE ====', flush=True)
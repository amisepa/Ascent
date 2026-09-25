# diagnose_specparam_diff.py — why do offsets disagree on synthetic spectra?
import os
import numpy as np
import scipy.io as sio

root = os.path.expandvars(r'C:\Users\ccann\Documents\MATLAB\Ascent')
M = sio.loadmat(os.path.join(root, 'validation_specparam', 'matlab_results.mat'), squeeze_me=True)
P = sio.loadmat(os.path.join(root, 'validation_specparam', 'python_results.mat'), squeeze_me=True)

npkM = np.asarray(M['npkM'], float).ravel()
npkP = np.asarray(P['npkPS'], float).ravel()
npkT = np.asarray(M['nPkTrue'], float).ravel()
doff = np.asarray(M['offMT'], float).ravel() - np.asarray(P['offPS'], float).ravel()
dexp = np.asarray(M['expMT'], float).ravel() - np.asarray(P['expPS'], float).ravel()

print('true peak count dist  :', np.bincount(npkT.astype(int), minlength=4))
print('ASCENT peak count dist:', np.bincount(npkM.astype(int), minlength=4))
print('specparam pk count    :', np.bincount(npkP.astype(int), minlength=4))

for t in range(4):
    m = npkT == t
    if m.sum():
        print(f'true peaks={t} (n={m.sum()}): dOff mean={doff[m].mean():+.3f} sd={doff[m].std():.3f} | dExp mean={dexp[m].mean():+.3f}')

# does offset error scale with how many peaks were MISSED?
missed = npkT - npkM
m = missed >= 0
print('corr(dOff, missed peaks) =', np.corrcoef(doff[m], missed[m])[0,1].round(3))
print('corr(dExp, missed peaks) =', np.corrcoef(dexp[m], missed[m])[0,1].round(3))

# spectra with NO true peaks: pure 1/f+noise. Compare offsets there (cleanest parity test)
z = npkT == 0
if z.sum():
    print(f'pure-1/f spectra (n={z.sum()}): dOff mean={doff[z].mean():+.4f} sd={doff[z].std():.4f} | dExp mean={dexp[z].mean():+.4f} sd={dexp[z].std():.4f}')
    print('  -> on peak-free spectra the two implementations should agree to noise level')

# real EEG peak counts
npkMR = np.asarray(M['npkM'], float).ravel()  # synthetic
print('real M npk mean:', np.nanmean(npkMR), '(synthetic col)')  # sanity only
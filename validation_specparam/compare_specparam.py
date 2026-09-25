# compare_specparam.py — agreement stats + supplementary figure for ASCENT vs specparam
import os, sys
import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root   = os.path.expandvars(r'C:\Users\ccann\Documents\MATLAB\Ascent')
outdir = os.path.join(root, 'validation_specparam')
mfile = sys.argv[1] if len(sys.argv) > 1 else 'matlab_results.mat'
M = sio.loadmat(os.path.join(root, 'validation_specparam', mfile), squeeze_me=True)
P = sio.loadmat(os.path.join(root, 'validation_specparam', 'python_results.mat'), squeeze_me=True)

def stats(a, b, mask=None):
    if mask is None:
        mask = np.isfinite(a) & np.isfinite(b)
    a, b = np.asarray(a)[mask], np.asarray(b)[mask]
    if len(a) < 2:
        return dict(n=len(a), r=np.nan, mae=np.nan, bias=np.nan, loa=np.nan)
    r = np.corrcoef(a, b)[0, 1]
    d = a - b
    bias, sd = d.mean(), d.std(ddof=1)
    return dict(n=len(a), r=r, mae=np.mean(np.abs(d)), bias=bias, loa=1.96*sd)

# ---- synthetic: exponent & offset ----
s_exp  = stats(M['expMT'], P['expPS'])
s_off  = stats(M['offMT'], P['offPS'])
s_r2   = stats(M['r2M'],   P['r2PS'])
s_mae  = stats(M['maeM'],  P['maePS'])
# peak-count agreement
npkM = np.asarray(M['npkM'], float).ravel(); npkP = np.asarray(P['npkPS'], float).ravel()
pk_same = np.mean(npkM == npkP) * 100
# ---- real: exponent & offset ----
r_exp = stats(M['expMR'], P['expPR'])
r_off = stats(M['offMR'], P['offPR'])

print('=== SYNTHETIC (n=200) ===')
print(f"exponent : r={s_exp['r']:.4f}  MAE={s_exp['mae']:.4f}  bias={s_exp['bias']:+.4f} +/- {s_exp['loa']:.4f} (95% LoA)")
print(f"offset   : r={s_off['r']:.4f}  MAE={s_off['mae']:.4f}  bias={s_off['bias']:+.4f} +/- {s_off['loa']:.4f}")
print(f"fit MAE  : r={s_mae['r']:.4f}   fit R2 : r={s_r2['r']:.4f}")
print(f"peak-count identical: {pk_same:.1f}%  (M mean {np.nanmean(npkM):.2f} vs P mean {np.nanmean(npkP):.2f})")
print('=== REAL (n=64) ===')
print(f"exponent : r={r_exp['r']:.4f}  MAE={r_exp['mae']:.4f}  bias={r_exp['bias']:+.4f} +/- {r_exp['loa']:.4f}")
print(f"offset   : r={r_off['r']:.4f}  MAE={r_off['mae']:.4f}  bias={r_off['bias']:+.4f} +/- {r_off['loa']:.4f}")

# ground-truth recovery (synthetic)
gt_exp = stats(np.asarray(M['expTrue'],float).ravel(), np.asarray(M['expMT'],float).ravel())
gt_off = stats(np.asarray(M['offTrue'],float).ravel(), np.asarray(M['offMT'],float).ravel())
gt_expP = stats(np.asarray(M['expTrue'],float).ravel(), np.asarray(P['expPS'],float).ravel())
print(f"ground-truth recovery  M: exp r={gt_exp['r']:.3f} MAE={gt_exp['mae']:.3f} | off r={gt_off['r']:.3f} MAE={gt_off['mae']:.3f}")
print(f"ground-truth recovery  P: exp r={gt_expP['r']:.3f} MAE={gt_expP['mae']:.3f}")

# ================= supplementary figure =================
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
plt.rcParams.update({'font.size': 9})

def scatter_panel(ax, a, b, title, lims=None):
    m = np.isfinite(a) & np.isfinite(b)
    ax.scatter(np.asarray(b)[m], np.asarray(a)[m], s=12, alpha=0.6, c='#3060c0', edgecolors='none')
    lim = lims or [np.nanmin([a[m], b[m]]), np.nanmax([a[m], b[m]])]
    ax.plot(lim, lim, 'k--', lw=1)
    r = np.corrcoef(np.asarray(a)[m], np.asarray(b)[m])[0, 1]
    ax.set_title(f'{title}\nr = {r:.3f}')
    ax.set_xlabel('specparam (Python)'); ax.set_ylabel('ASCENT (MATLAB)')
    ax.set_xlim(lim); ax.set_ylim(lim)

def ba_panel(ax, a, b, title):
    m = np.isfinite(a) & np.isfinite(b)
    a, b = np.asarray(a)[m], np.asarray(b)[m]
    mean_ab, d = (a+b)/2, a-b
    bias, sd = d.mean(), d.std(ddof=1)
    ax.scatter(mean_ab, d, s=12, alpha=0.6, c='#c06018', edgecolors='none')
    ax.axhline(bias, color='k', lw=1)
    ax.axhline(bias+1.96*sd, color='gray', ls='--', lw=1)
    ax.axhline(bias-1.96*sd, color='gray', ls='--', lw=1)
    ax.set_title(f'{title}\nbias {bias:+.3f}, LoA ±{1.96*sd:.3f}')
    ax.set_xlabel('mean of both'); ax.set_ylabel('ASCENT − specparam')

scatter_panel(axes[0,0], M['expMT'], P['expPS'], 'Synthetic exponent')
scatter_panel(axes[0,1], M['offMT'], P['offPS'], 'Synthetic offset')
ba_panel(axes[0,2], M['expMT'], P['expPS'], 'Synthetic exponent (B-A)')
ba_panel(axes[0,3], M['offMT'], P['offPS'], 'Synthetic offset (B-A)')
scatter_panel(axes[1,0], M['expMR'], P['expPR'], 'Real EEG exponent (64 ch)')
scatter_panel(axes[1,1], M['offMR'], P['offPR'], 'Real EEG offset (64 ch)')
ba_panel(axes[1,2], M['expMR'], P['expPR'], 'Real exponent (B-A)')
ba_panel(axes[1,3], M['offMR'], P['offPR'], 'Real offset (B-A)')
fig.suptitle('ASCENT (MATLAB) vs specparam 2.0.0rc6 (Python) — matched settings, fixed mode, 1–40 Hz', fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
png = os.path.join(outdir, sys.argv[2] if len(sys.argv) > 2 else 'specparam_validation_figure.png')
fig.savefig(png, dpi=300)
print('saved', png)
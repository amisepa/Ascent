# Lessons

- MATLAB graphics (figure + drawnow/print) hang under `matlab -batch`, `-r` and
  `-noFigureWindows` when launched from the agent shell on this machine; each
  attempt leaves MATLABWindow.exe processes. Numeric runs work headless. Make
  figures from the MATLAB desktop (make_revision_figures.m).
- In -batch runs, a second `parfor` call (eeg_robust_statistics) can hang;
  use serial copies of run_stats_permutation/compute_mcc for headless reruns.
- specparam port: a median+2*MAD robust initial fit keeps large peak humps and
  pulls the aperiodic line up, so peaks go undetected and offsets inflate. Use
  specparam's fit (OLS, clip residuals at 0, refit on points at/below the line).
- ascent_compute auto-saves the .set to EEG.filepath: redirect EEG.filepath to
  a scratch folder when batch-processing data you must not modify.
- Word COM edits under Track Changes: build equations in a scratch paragraph at
  the end and paste (BuildUp mid-text eats characters); turn off
  Options.SmartCutPaste (it deletes spaces next to replaced equations). Google
  Docs can't import tracked edits inside equations or Word's "#(n)" numbering.

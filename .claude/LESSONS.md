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
- Manuscript figures: a positive w:ind firstLine (MDPI 2608) on a figure
  paragraph pushes a 5-in image off the page; set it to 0 but keep explicit
  firstLine="0" (it overrides the body style's first-line indent). The v13 PDF
  was rendered at ~74% and hid this: check image bboxes in the new PDF (pymupdf).
- plot_clusters(summary, mask, tvals, tvals, ...) plots t ± 1.96 under whatever
  label is passed; give it the per-subject difference [chan x freq x subj] to
  get mean ± 95% CI in data units.
- docx_xml.replace_untracked gives the whole replacement the first run's
  formatting: keep the span inside one run ("Fig 10." not "Fig 10. Differences").
- compute_mcc forms clusters on |t| and pull_clusters (2D path) labels the
  unsigned mask, so touching positive/negative regions become one cluster named
  after its peak sign (MFE/RCMFE mean lost their fine-scale negative cluster).
  Report with pull_clusters_by_sign (make_revision_figures.m,
  ascent_group_analysis.m) and seed each test with rng(1).
- The time-resolved aperiodic fit in make_revision_figures runs ~10 min per
  domain; it only looked hung because 'progress' was off (now on, plus
  timestamps). Other MATLAB batch jobs on this machine (HBN ICA, 12 cores) slow it.

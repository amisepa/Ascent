%% Save High-res GUI figures for tutorial
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).


%% Normal figures

cd 'C:\Users\ccann\Documents\MATLAB\Ascent\figures_new2\subject_level'
saveas(gcf, 'fig_rcmfe_ica_dipole.fig');
print(gcf, 'fig_rcmfe_ica_dipole.png', '-dpng', '-r300');


%% GUI figures (for Matlab 2025 or later)

cd 'C:\Users\ccann\Documents\MATLAB\Ascent\figures_new2'
saveas(gcf, 'fig1c.fig');
exportapp(gcf,'fig1c.png')


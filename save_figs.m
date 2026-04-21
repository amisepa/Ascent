%% Save High-res GUI figures for tutorial


%% Normal figures

cd 'C:\Users\ccann\Documents\MATLAB\Ascent\figures_new\subject_level'
saveas(gcf, 'fig_rcmfe_mean_ica.fig');
print(gcf, 'fig_rcmfe_mean_ica.png', '-dpng', '-r300');



%% GUI figures (for Matlab 2025 or later)

cd 'C:\Users\ccann\Documents\MATLAB\Ascent\figures_new\subject_level'
saveas(gcf, 'fig1c.fig');
exportapp(gcf,'fig1c.png')


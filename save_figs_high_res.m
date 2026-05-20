%% Save High-res GUI figures for tutorial


%% Normal figures

cd 'C:\Users\ccann\Documents\MATLAB\Ascent\figures_new2\subject_level'
saveas(gcf, 'fig_rcmfe_std.fig');
print(gcf, 'fig_rcmfe_std.png', '-dpng', '-r300');



%% GUI figures (for Matlab 2025 or later)

cd 'C:\Users\ccann\Documents\MATLAB\Ascent\figures_new'
saveas(gcf, 'fig1a.fig');
exportapp(gcf,'fig1a.png')


function ascent_plot_dipoles(vals, icinfo, dipfit, measure, fs, scales, rv_thresh)
% ascent_plot_dipoles  Equivalent dipoles of the ICs, colored by a measure, and the
% measure by brain lobe.
%
% Called by ascent_plot_multiscale and ascent_plot_aperiodic in ICA mode when the
% dataset has DIPFIT models; not intended for direct use.
%   vals     [ICs x scales] (multiscale measures) or [ICs x 1] (one value per IC)
%   icinfo   per-IC table from ascent_ic_info (ICLabel class, RV, position, AAL lobe)
%   dipfit   EEG.dipfit
%   measure  name for the titles and color bar (e.g. 'RCMFE', 'Aperiodic exponent')
%   fs       font size (points)
%   scales   scale values of the columns of vals (multiscale only)
%
% Dipolar ICs (RV < rv_thresh, default 0.15) located inside the brain are drawn on
% the MNI template: color = the measure (its mean across scales for multiscale
% measures), size = fit quality (1 - RV), circle = ICLabel brain, square = other
% ICLabel classes. Dipoles outside the brain are artifact ICs (eyes, muscle) and are
% left out. Below, the measure by AAL lobe for lobes with at least 3 ICs (with
% scale-resolved curves for multiscale measures). A table of the dipoles is printed
% in the command window.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).

if nargin < 6, scales = []; end
if nargin < 7 || isempty(rv_thresh), rv_thresh = 0.15; end
if isempty(icinfo) || isempty(dipfit), return; end
if ~isfield(dipfit, 'coordformat') || ~strcmpi(dipfit.coordformat, 'MNI'), return; end   % warned by ascent_ic_info
lw = max(1, fs/9);
n_ics = size(vals, 1);
multi = size(vals, 2) > 1;
if multi
    if isempty(scales) || ~isnumeric(scales), scales = 1:size(vals, 2); end
    scales = scales(:)';
    cName = ['mean ' measure];
else
    cName = measure;
end
mean_all = mean(vals, 2, 'omitnan');
sel = icinfo([icinfo.dipolar] & [icinfo.inbrain] & [icinfo.ic] <= n_ics);
if ~isempty(sel), sel = sel(isfinite(mean_all([sel.ic]))); end
if isempty(sel)
    warning('ascent_plot: no dipolar IC (RV < %.0f%%) inside the brain; no dipole figure.', rv_thresh*100);
    return
end
val_sel   = vals([sel.ic],:);
mean_v    = mean_all([sel.ic]);
hasICL    = ~all(cellfun(@isempty, {sel.iclabel}));
isBrainIC = strcmp({sel.iclabel}, 'Brain');

% Table in the command window
fprintf('%s, dipolar ICs inside the brain (RV < %.0f%%): %d\n', measure, rv_thresh*100, numel(sel));
fprintf('   IC     RV   MNI x, y, z (mm)   AAL region (lobe)                      ICLabel (p)          %s\n', ...
        [upper(cName(1)) cName(2:end)]);
for i = 1:numel(sel)
    s = sel(i);
    if isempty(s.aal), reg = 'unlabeled'; else, reg = sprintf('%s (%s)', s.aal, s.lobe); end
    if isempty(s.iclabel), icl = '-'; else, icl = sprintf('%s (%.2f)', s.iclabel, s.p_iclabel); end
    fprintf('  %3d  %4.1f%%  %5.0f %5.0f %5.0f    %-38s %-20s %.3f\n', s.ic, 100*s.rv, s.mni, reg, icl, mean_v(i));
end

% Lobes with at least min_n ICs
lobe_names = {'Frontal','Central','Parietal','Temporal','Occipital','Cingulate','Insula','Subcortical','Cerebellum'};
lobe_clrs  = {[0.00 0.45 0.74],[0.90 0.55 0.00],[0.85 0.20 0.55],[0.45 0.30 0.85],[0.45 0.45 0.45], ...
              [0.15 0.15 0.55],[0.30 0.75 0.93],[0.00 0.00 0.00],[0.60 0.40 0.20]};
min_n = 3;
ic_lobes = {sel.lobe};
has_reg  = any(cellfun(@(r) sum(strcmp(ic_lobes, r)) >= min_n, lobe_names));

% ---- Run dipplot into temp figure, capture as image --------------------
dip_img = [];
v_row = mean_v';
cmap = plasma_map(256);   % purple to yellow, no green, readable on the grey MRI
cl = [min(v_row) max(v_row)];
if ~(cl(2) > cl(1)), cl = cl(1) + [-1 1] * max(1e-3 * abs(cl(1)), 1e-6); end   % one value (or all equal)
if isfield(dipfit,'mrifile')
    try
        % Cluster membership is not marked: few high-entropy ICs are dipolar,
        % and a dipole position only means something for those that are.
        sources = struct('posxyz', {sel.mni}, 'momxyz', [0 0 0], 'rv', {sel.rv});
        ci = 1 + round((v_row - cl(1)) / (cl(2) - cl(1)) * 255);
        colors = arrayfun(@(i) cmap(i,:), ci, 'UniformOutput', false);
        msz = 5 + 7 * (1 - min([sel.rv], rv_thresh) / rv_thresh);   % RV 0: 12 pt; RV at threshold: 5 pt
        mrk = repmat({'o'}, 1, numel(sel));
        if hasICL, mrk(~isBrainIC) = {'s'}; end
        pre_figs = findall(0,'Type','figure');
        dipplot(sources,'mri',template_mri_without_slab(dipfit),'coordformat','MNI','color',colors,...
            'summary','3d','dipolesize',15,'dipolelength',0,'gui','off','mesh','off');
        drawnow; pause(0.3);
        new_figs = setdiff(findall(0,'Type','figure'),pre_figs);
        if ~isempty(new_figs)
            tmp_fig = new_figs(end);
            hide_edge_on_slices(tmp_fig);
            outline_dipoles(tmp_fig, msz, mrk);
            set(tmp_fig,'Color','k','InvertHardcopy','off');
            try
                dip_img = print(tmp_fig,'-RGBImage','-r300');  % sharper than getframe
            catch
                dip_img = getframe(tmp_fig).cdata;
            end
            close(tmp_fig);
        end
    catch ME
        warning('ascent_plot: dipplot capture failed - %s', ME.message);
    end
end
has_dip = ~isempty(dip_img);
if ~has_dip && ~has_reg, return; end

% ---- Combined figure layout --------------------------------------------
fig = figure('Color','w','InvertHardCopy','off','NumberTitle','off',...
             'Name',['Dipole sources + ' measure ' by region']);
try icadefs; set(fig,'color',BACKCOLOR); catch; end %#ok<NODEF>
set(fig,'Units','pixels'); p0 = get(fig,'Position');
if has_dip && has_reg
    sz = [1000 720]; rowH = 0.50; dipPos = [0.01 0.52 0.88 0.38];
elseif has_dip
    sz = [1000 380]; dipPos = [0.01 0.03 0.88 0.72];
else
    sz = [1000 420]; rowH = 1;
end
set(fig,'Position',[p0(1) max(40, p0(2)+p0(4)-sz(2)) sz]);
if has_dip
    ax_dip = axes(fig,'Position',dipPos);   % room on the right for the color bar and its label
    image(ax_dip, dip_img); axis(ax_dip,'image','off');
    colormap(ax_dip, cmap); clim(ax_dip, cl);
    cb = colorbar(ax_dip, 'eastoutside'); cb.FontSize = fs; cb.FontWeight = 'bold';
    cb.Label.String = [upper(cName(1)) cName(2:end)]; cb.Label.FontSize = fs; cb.Label.FontWeight = 'bold';
    cb.Label.Interpreter = 'none';
    t2 = 'Size: fit quality (1 - RV)';
    if hasICL
        t2 = sprintf('%s; circles: ICLabel brain (%d), squares: other classes (%d)', t2, ...
                     sum(isBrainIC), sum(~isBrainIC));
    end
    title(ax_dip, {sprintf('Equivalent dipoles inside the brain (RV < %.0f%%, %d ICs), color: %s', ...
          rv_thresh*100, numel(sel), cName), t2}, 'FontSize',fs,'FontWeight','bold','Interpreter','none');
end
if ~has_reg
    fprintf('No AAL lobe has %d or more dipolar ICs: no region plot.\n', min_n);
    return
end
% OuterPosition: MATLAB fits the tick labels, axis labels and titles inside
if multi
    ax1 = axes(fig,'OuterPosition',[0    0 0.64 rowH]);
    ax2 = axes(fig,'OuterPosition',[0.64 0 0.36 rowH]);
    hold(ax1,'on'); box(ax1,'on');
else
    ax1 = gobjects(0);
    ax2 = axes(fig,'OuterPosition',[0.25 0 0.50 rowH]);
end
hold(ax2,'on'); box(ax2,'on');

xpos=0; xt_pos=[]; xt_lbl={}; leg_h=gobjects(0); allv=[];
for rr = 1:numel(lobe_names)
    rname = lobe_names{rr};
    idx   = strcmp(ic_lobes, rname);
    n_r   = sum(idx);
    if n_r < min_n, continue; end
    clr = lobe_clrs{rr};

    if multi
        mu  = mean(val_sel(idx,:),1,'omitnan');
        sem = std( val_sel(idx,:),0,1,'omitnan') / sqrt(n_r);
        fill(ax1,[scales fliplr(scales)],[mu+sem fliplr(mu-sem)],...
             clr,'FaceAlpha',0.15,'EdgeColor','none');
        leg_h(end+1) = plot(ax1,scales,mu,'Color',clr,'LineWidth',2.2*lw,...
            'DisplayName',sprintf('%s (n=%d)',rname,n_r)); %#ok<AGROW>
    end

    xpos = xpos+1;
    v = mean_v(idx); allv = [allv; v]; %#ok<AGROW>
    scatter(ax2,xpos+(rand(n_r,1)-0.5)*0.3,v,30*lw^2,clr,...
            'filled','MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
    plot(ax2,xpos+[-0.3 0.3],repmat(mean(v,'omitnan'),1,2),...
         '-','Color',clr*0.6,'LineWidth',2.5*lw);
    if multi
        xt_lbl{end+1} = rname;                         %#ok<AGROW> (n is in the legend)
    else
        xt_lbl{end+1} = sprintf('%s (n=%d)', rname, n_r); %#ok<AGROW>
    end
    xt_pos(end+1) = xpos; %#ok<AGROW>
end

if multi
    xlim(ax1,[scales(1) scales(end)]);
    xlabel(ax1,'Scale'); ylabel(ax1,measure,'Interpreter','none');
    title(ax1,sprintf('Scale-resolved %s by AAL lobe  (mean %s SEM)', measure, char(177)),'Interpreter','none');
    set(ax1,'TickDir','out');
end
xlim(ax2,[0 xpos+1]);
allv = allv(isfinite(allv));
if numel(allv) > 1 && max(allv) > min(allv)   % pad so no point sits on the frame
    pad = 0.08*(max(allv)-min(allv)); ylim(ax2,[min(allv)-pad max(allv)+pad]);
end
set(ax2,'XTick',xt_pos,'XTickLabel',xt_lbl,'TickDir','out','XTickLabelRotation',30);
ylabel(ax2,[upper(cName(1)) cName(2:end)],'Interpreter','none');
title(ax2,[upper(cName(1)) cName(2:end) ' by AAL lobe'],'Interpreter','none');
set([ax1 ax2],'FontSize',fs,'FontWeight','bold', ...
    'TitleFontSizeMultiplier',1,'LabelFontSizeMultiplier',1);
if multi, legend(ax1,leg_h,'Location','best','Box','off','FontSize',max(fs-1,8)); end
end


%% Purple-to-yellow colour map (anchors of matplotlib's plasma)
function c = plasma_map(n)
a = [0.050 0.030 0.528; 0.417 0.000 0.658; 0.692 0.165 0.565; 0.881 0.393 0.383; ...
     0.988 0.652 0.211; 0.940 0.975 0.131];
c = interp1(linspace(0,1,size(a,1)), a, linspace(0,1,n));
end


%% Dipole markers with a black edge (dots vanish on the grey MRI)
function outline_dipoles(fig, msz, mrk)
% Size and shape per dipole; dipplot tags each marker 'dipole<index>'
for h = findall(fig, 'Type', 'line', 'Marker', '.')'
    t = get(h,'Tag');
    if ~strncmp(t, 'dipole', 6), continue; end
    i = str2double(t(7:end));
    if isnan(i) || i > numel(msz), continue; end
    c = get(h, 'Color');
    set(h, 'Marker', mrk{i}, 'MarkerSize', msz(i) * ifelse(strcmp(mrk{i},'s'), 0.9, 1), ...
        'MarkerFaceColor', c, 'MarkerEdgeColor', [0 0 0], 'LineWidth', 1);
end
end


%% dipplot '3d' summary: show only the slice each view faces
function hide_edge_on_slices(fig)
% Each of the three views shows one MRI slice face-on; hide the two seen
% edge-on and use exact axis views with an orthographic projection. (The bar
% above the coronal view was not one of these: see template_mri_without_slab.)
for ax = findall(fig,'Type','axes')'
    hs = findall(ax,'Type','surface','Tag','img');
    if isempty(hs), continue; end
    v = get(ax,'CameraPosition') - get(ax,'CameraTarget');
    [~, d] = max(abs(v));                        % coordinate the view looks along
    e = zeros(1,3); e(d) = sign(v(d)); view(ax, e);
    set(ax,'Projection','orthographic');
    w = zeros(numel(hs),1);
    for i = 1:numel(hs)
        c = get(hs(i), {'XData','YData','ZData'}); c = c{d};
        w(i) = max(c(:)) - min(c(:));             % slice extent along the view
    end
    [~, keep] = min(w);
    set(hs(setdiff(1:numel(hs), keep)), 'Visible','off');
end
end


%% EEGLAB template MRI without the stray tissue above the scalp
function mri = template_mri_without_slab(dipfit)
% The top slices of EEGLAB's template MRI (dipfit standard_mri.mat) hold stray
% tissue 5-8 mm above the scalp, which dipplot draws as a bar over the coronal
% view and blobs over the sagittal one. Blank every voxel above 60 mm that lies
% more than 4 mm above the scalp of the BEM head model. Any other case (other
% MRI or head model, rotated voxel grid, read error) returns the file name, so
% dipplot shows the MRI as before.
mri = dipfit.mrifile;
if exist(mri, 'file') ~= 2                    % stored without its folder: the dipfit plugin's file
    [~, n, e] = fileparts(mri); c = fullfile(fileparts(which('dipfitdefs')), 'standard_BEM', [n e]);
    if exist(c, 'file') == 2, mri = c; end
end
try
    M = load('-mat', mri);  A = M.mri.anatomy;  T = M.mri.transform;
    if any(any(abs(T(1:3,1:3) - diag(diag(T(1:3,1:3)))) > 1e-9)), return; end
    hdm = dipfit.hdmfile;
    if exist(hdm, 'file') ~= 2
        [~, n, e] = fileparts(hdm); hdm = fullfile(fileparts(which('dipfitdefs')), 'standard_BEM', [n e]);
    end
    H = load('-mat', hdm);  P = [];
    for b = H.vol.bnd(:)'                      % scalp: the outermost surface
        if isfield(b,'pos'), q = b.pos; else, q = b.pnt; end
        if isempty(P) || max(q(:,3)) > max(P(:,3)), P = q; end
    end
    x = T(1,1)*(1:size(A,1)) + T(1,4);         % voxel centres in mm (1-based)
    y = T(2,2)*(1:size(A,2)) + T(2,4);
    z = T(3,3)*(1:size(A,3)) + T(3,4);
    Q = sortrows(P(P(:,3) > 0, :), -3);       % upper half; highest point first
    [~, iu] = unique(Q(:,1:2), 'rows', 'stable');  Q = Q(iu,:);   % one point per column
    [X, Y] = ndgrid(x, y);
    ztop = griddata(Q(:,1), Q(:,2), Q(:,3), X, Y, 'linear');   % scalp height per column
    ztop(isnan(ztop)) = -Inf;                  % outside the scalp footprint
    A(bsxfun(@gt, reshape(z,1,1,[]), max(ztop + 4, 60))) = 0;
    mri = M.mri;  mri.anatomy = A;
catch
end
end


function out = ifelse(cond,a,b)
if cond, out = a; else, out = b; end
end

function info = ascent_ic_info(nIC, iclabel, dipfit, chanlocs, rv_thresh)
% ascent_ic_info - ICLabel class and equivalent dipole of each IC, in one table.
%
%   info = ascent_ic_info(size(EEG.icaweights,1), EEG.etc.ic_classification.ICLabel, EEG.dipfit)
%   info = ascent_ic_info(nIC, iclabel, dipfit, EEG.chanlocs(EEG.icachansind), 0.15)
%
% iclabel is EEG.etc.ic_classification.ICLabel (or its [nIC x 7] classifications
% matrix, classes in ICLabel's order); dipfit is EEG.dipfit. Either may be empty.
% Returns a struct array, info(k) for IC k:
%   ic        - IC index
%   iclabel   - most probable ICLabel class ('' without ICLabel)
%   p_iclabel - probability of that class
%   p_brain   - probability of the Brain class
%   rv        - residual variance of the dipole fit (0-1; NaN without a fit):
%               the share of the IC map one dipole does not explain. It
%               measures how dipolar the map is, not how accurate the position is.
%   dipolar   - rv < rv_thresh (default 0.15, the EEGLAB convention)
%   mni       - dipole position (mm, MNI)
%   inbrain   - inside the innermost surface of the head model (BEM)
%   aal       - AAL region containing the dipole; when the position is not
%               labeled (e.g., white matter), the nearest labeled region within
%               10 mm, prefixed by 'near '; '' otherwise
%   lobe      - Frontal, Central, Parietal, Temporal (incl. hippocampus,
%               amygdala), Occipital, Cingulate, Insula, Subcortical or
%               Cerebellum (from the AAL region); '' when unlabeled
%
% Labels use the AAL atlas (Tzourio-Mazoyer et al., 2002) shipped with
% FieldTrip in EEGLAB, by containment rather than by distance to a centroid
% (centroid rules mislabel midline and deep dipoles). With chanlocs, it also
% warns when DIPFIT used an identity coord_transform with channel locations in
% EEGLAB's convention (nose along +X): the montage is then rotated by 90 deg
% relative to the MNI head model and the dipoles must be refit (for positions
% from standard_1005: coord_transform [0 0 0 0 0 -pi/2 1 1 1]).
%
% ICLabel was trained on average-referenced data decomposed with Infomax;
% treat its classes as a guide when the pipeline differs.
%
% Copyright (C) 2021-2026 Cedric Cannard. ASCENT EEGLAB plugin, GNU GPL v3 or later (see LICENSE).

if nargin < 2, iclabel = []; end
if nargin < 3, dipfit = []; end
if nargin < 4, chanlocs = []; end
if nargin < 5 || isempty(rv_thresh), rv_thresh = 0.15; end

% ICLabel
P = []; classes = {'Brain','Muscle','Eye','Heart','Line Noise','Channel Noise','Other'};
if isstruct(iclabel) && isfield(iclabel, 'classifications')
    P = iclabel.classifications;
    if isfield(iclabel, 'classes') && ~isempty(iclabel.classes), classes = iclabel.classes; end
elseif isnumeric(iclabel)
    P = iclabel;
end
if ~isempty(P) && size(P,1) ~= nIC
    warning('ascent: ICLabel has %d ICs, the data %d; ICLabel ignored.', size(P,1), nIC); P = [];
end
iBrain = find(strcmpi(classes, 'Brain'), 1);

hasDip = isstruct(dipfit) && isfield(dipfit, 'model') && ~isempty(dipfit.model);
if hasDip
    if ~isempty(chanlocs), check_frame(dipfit, chanlocs); end
    inBrain = brain_test(dipfit);
    [A, names] = load_aal();
end

info = struct('ic', num2cell(1:nIC), 'iclabel', '', 'p_iclabel', NaN, 'p_brain', NaN, ...
    'rv', NaN, 'dipolar', false, 'mni', [NaN NaN NaN], 'inbrain', false, 'aal', '', 'lobe', '');
for k = 1:nIC
    if ~isempty(P)
        [pk, c] = max(P(k,:));
        info(k).iclabel = classes{c}; info(k).p_iclabel = pk;
        if ~isempty(iBrain), info(k).p_brain = P(k, iBrain); end
    end
    if ~hasDip || k > numel(dipfit.model), continue; end
    m = dipfit.model(k);
    if isempty(m.posxyz) || isempty(m.rv), continue; end
    p = m.posxyz(1,:);
    [info(k).aal, info(k).lobe] = aal_label(p, A, names);
    info(k).rv = m.rv; info(k).dipolar = m.rv < rv_thresh;
    info(k).mni = p; info(k).inbrain = inBrain(p);
end
end

%% Inside the innermost BEM surface
function f = brain_test(dipfit)
f = @(p) true;
try
    H = load('-mat', dipfit.hdmfile);
    B = H.vol.bnd(:); ext = zeros(numel(B),1); P = cell(numel(B),1);
    for i = 1:numel(B)
        if isfield(B(i),'pos'), P{i} = B(i).pos; else, P{i} = B(i).pnt; end
        ext(i) = max(max(P{i}) - min(P{i}));
    end
    [~, k] = min(ext);
    DT = delaunayTriangulation(P{k});
    f = @(p) ~isnan(pointLocation(DT, p(:)'));
catch
end
end

%% AAL atlas (NIfTI-1, big-endian int16) read without toolboxes, cached
function [A, names] = load_aal()
persistent cache
if ~isempty(cache), A = cache.A; names = cache.names; return; end
A = []; names = containers.Map('KeyType','double','ValueType','char');
f = '';
ft = which('ft_defaults');
if ~isempty(ft), f = fullfile(fileparts(ft), 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii'); end
if isempty(f) || ~exist(f, 'file')   % FieldTrip-lite not on the path: look in EEGLAB's plugins
    d = dir(fullfile(fileparts(which('eeglab')), 'plugins', 'Fieldtrip*', 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii'));
    if ~isempty(d), f = fullfile(d(1).folder, d(1).name); end
end
if isempty(f) || ~exist(f, 'file'), return; end
try
    fid = fopen(f, 'r', 'ieee-be');
    if fread(fid, 1, 'int32') ~= 348, fclose(fid); fid = fopen(f, 'r', 'ieee-le'); fread(fid, 1, 'int32'); end
    fseek(fid, 40, 'bof'); dim = fread(fid, 8, 'int16')';
    fseek(fid, 70, 'bof'); dtype = fread(fid, 1, 'int16');
    fseek(fid, 108, 'bof'); off = fread(fid, 1, 'float32');
    fseek(fid, 280, 'bof'); S = reshape(fread(fid, 12, 'float32'), 4, 3)';
    if dtype ~= 4, fclose(fid); return; end            % int16 expected
    fseek(fid, round(off), 'bof');
    V = reshape(fread(fid, prod(dim(2:4)), 'int16=>double'), dim(2:4));
    fclose(fid);
    A = struct('V', V, 'S', [S; 0 0 0 1]);
    T = readcell(strrep(f, '.nii', '.txt'), 'Delimiter', '\t', 'FileType', 'text');
    for i = 1:size(T,1), names(double(T{i,3})) = T{i,2}; end
    cache = struct('A', A, 'names', names);
catch
    A = [];
end
end

function [aal, lobe] = aal_label(p, A, names)
aal = ''; lobe = '';
if isempty(A), return; end
ijk = A.S \ [p(:); 1];                       % 0-based voxel coordinates
v = round(ijk(1:3))' + 1;
sz = size(A.V);
code = 0;
if all(v >= 1 & v <= sz), code = A.V(v(1), v(2), v(3)); end
prefix = '';
if code == 0                                  % not labeled: nearest labeled voxel within 10 mm (5 voxels)
    r = 5; [di, dj, dk] = ndgrid(-r:r, -r:r, -r:r);
    pts = [di(:) dj(:) dk(:)] + v; d = sqrt(sum(([di(:) dj(:) dk(:)] * 2).^2, 2));
    ok = all(pts >= 1 & pts <= sz, 2) & d <= 10;
    pts = pts(ok,:); d = d(ok);
    c = A.V(sub2ind(sz, pts(:,1), pts(:,2), pts(:,3)));
    if any(c > 0)
        [~, i] = min(d + 1e6 * (c == 0)); code = c(i); prefix = 'near ';
    end
end
if code > 0 && isKey(names, code)
    aal = [prefix names(code)];
    lobe = aal_lobe(names(code));
end
end

function lobe = aal_lobe(n)
if startsWith(n, {'Precentral', 'Postcentral', 'Paracentral', 'Rolandic'}), lobe = 'Central';
elseif startsWith(n, {'Frontal', 'Supp_Motor', 'Olfactory', 'Rectus'}), lobe = 'Frontal';
elseif startsWith(n, {'Parietal', 'SupraMarginal', 'Angular', 'Precuneus'}), lobe = 'Parietal';
elseif startsWith(n, {'Temporal', 'Heschl', 'Fusiform', 'Hippocampus', 'ParaHippocampal', 'Amygdala'}), lobe = 'Temporal';
elseif startsWith(n, {'Occipital', 'Calcarine', 'Cuneus', 'Lingual'}), lobe = 'Occipital';
elseif startsWith(n, 'Cingulum'), lobe = 'Cingulate';
elseif startsWith(n, 'Insula'), lobe = 'Insula';
elseif startsWith(n, {'Caudate', 'Putamen', 'Pallidum', 'Thalamus'}), lobe = 'Subcortical';
elseif startsWith(n, {'Cerebel', 'Vermis'}), lobe = 'Cerebellum';
else, lobe = 'Other';
end
end

%% Warn when the montage is rotated relative to the MNI head model
function check_frame(dipfit, chanlocs)
try
    if ~isfield(dipfit, 'coordformat') || ~strcmpi(dipfit.coordformat, 'MNI'), return; end
    t = dipfit.coord_transform;
    if isempty(t) || numel(t) < 6 || any(abs(t(4:6)) > 1e-6), return; end   % some rotation set
    lab = lower({chanlocs.labels});
    i = find(ismember(lab, {'fpz', 'fz', 'nz', 'afz'}), 1);
    if isempty(i), return; end
    if abs(chanlocs(i).X) > 2 * abs(chanlocs(i).Y) && chanlocs(i).X > 0
        warning(['ascent: DIPFIT used an identity coord_transform, but the channel locations ' ...
            'have the nose along +X (EEGLAB convention) while the MNI head model has it along +Y. ' ...
            'The dipoles were fit with the montage rotated by 90 deg: refit them (for positions ' ...
            'from standard_1005, coord_transform [0 0 0 0 0 -pi/2 1 1 1]).']);
    end
catch
end
end

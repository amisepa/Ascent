% eegplugin_ascent()
%
% ASCENT: Aperiodic Signal Complexity Estimation for Neurophysiological Time series
%
% Copyright (C) - EEGLAB, Cedric Cannard, August 2022-2025
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function vers = eegplugin_ascent(fig,try_strings,catch_strings)

% Plugin version
vers = '2.1';
if nargin < 3, return; end

% Add paths to subfolders
p = fileparts(which('eegplugin_ascent.m'));
addpath(p);
addpath(fullfile(p,'functions'))

% One 'ASCENT' entry in the EEGLAB menu bar that opens the ASCENT window directly
% (as BrainBeats does). Enabled at startup too: with no dataset loaded the window
% loads one, or the sample dataset (ascent_compute_gui). Epoched data work as
% well (ascent_compute concatenates epochs). The command goes to the history and
% the dataset with its outputs becomes a new EEGLAB dataset.
cb_compute = [ try_strings.no_check '[EEG,LASTCOM] = ascent_compute(EEG);' catch_strings.new_and_hist ];
uimenu(fig, 'Label', 'ASCENT', 'Tag', 'menu_ascent', ...
    'userdata', 'startup:on;continuous:on;epoch:on;study:off;erpset:off', ...
    'Callback', cb_compute);

end




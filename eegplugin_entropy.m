% eegplugin_heartbrain() - EEGLAB plugin for different types of entropy measures:
% sample entropy, multiscale entropy, refined composite multiscale entropy. 
% Implements the option to bandpass-filter each scale to remove spectral biases.
% 
% input: 
% 
% output: 
% 
% Usage: eegplugin_entropy(fig, trystrs, catchstrs);
%
% Inputs:
%   fig        - [integer]  EEGLAB figure
%   trystrs    - [struct] "try" strings for menu callbacks.
%   catchstrs  - [struct] "catch" strings for menu callbacks.
%
% Copyright (C) - Cedric Cannard, August 2022
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

function vers = eegplugin_entropy(fig, trystrs, catchstrs)

%plugin version
vers = 'pop_entropy1.0';

if nargin < 3
    error('eegplugin_entropy requires 3 arguments');
end

%Add folder to path
p = which('eegplugin_entropy.m');
p = p(1:strfind(p,'eegplugin_entropy.m')-1);
if ~exist('eegplugin_entropy','dir')
    addpath(p);
end

%Find menu to import data 
menui = findobj(fig, 'tag', 'import data');

%Menu callbacks
comcnt = [trystrs.no_check '[EEG, LASTCOM] = pop_entropy;'  catchstrs.new_non_empty];

%Create menus
uimenu(menui, 'label', 'file containing EEG data', 'separator', 'on', 'callback', comcnt);

end

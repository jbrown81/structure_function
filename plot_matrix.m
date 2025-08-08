function plot_matrix(M, ci, clim)
% PLOT_MATRIX  Visualise a connectivity / weight matrix with module stripes
%
%   plot_matrix(M)                         % clim defaults to [-1 1]
%   plot_matrix(M, ci)                     % supply module IDs
%   plot_matrix(M, ci, [lo hi])            % custom colour limits
%
%   Inputs
%   ------      (• = required, ° = optional)
%   • M    : n×n matrix to display
%   ° ci   : n×1 vector of module IDs. If omitted, the function looks for
%            a variable named 'ci_15mods_v2' in the caller's workspace.
%   ° clim : 1×2 vector [lo hi] giving colour limits (default = [-1 1]).

% Copyright (c) 2025 Radiata, Inc.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

% ----------------------- argument handling --------------------------------
if nargin < 2 || isempty(ci)
    if evalin('base','exist(''ci_15mods_v2'',''var'')')
        ci = evalin('base','ci_15mods_v2');
    else
        error('Module ID vector ''ci'' must be supplied.');
    end
end
if nargin < 3 || isempty(clim)
    clim = [-1 1];
end
n = size(M,1);

% ---------------------------- constants -----------------------------------
mod_names = {'VIS','SM','DA','AUD','HIP','PHC','AMY', ...
             'CO','SAL','AT','FPl','FPr','DMN','SUB'};

cols_14mods = [117  31 133;
                74 131 178;
               111 192 103;
               237 252 172;
               152 227 168;
                65 202 163;
                10 216 216;
               251 154 212;
               213  85 252;
               220 115 152;
               228 147  51;
               132  51  73;
               201  66  80;
                 0   0   0] ./ 255;

% simple blue‑white‑red map for the matrix
neg = [linspace(0,1,128)', linspace(0,1,128)', ones(128,1)];    % blue→white
pos = [ones(128,1), linspace(1,0,128)', linspace(1,0,128)'];    % white→red
bwr = [neg ; pos(2:end,:)];

% --------------------------- re‑ordering ----------------------------------
[~, ord] = sort(ci);
ci_sorted = ci(ord);

% --------------------------- figure setup ---------------------------------
figure('Position',[100 100 1200 1200]);

% ---- subplot (1): module colour bar on the left --------------------------
h1 = subplot(2,2,1);
imagesc(sort(ci));                     % a 1‑D colour stripe
set(h1,'XTick',[], ...
        'YTick',module_centres(ci_sorted), ...
        'YTickLabel',mod_names, ...
        'FontSize',18,'TickLength',[0 0]);
h1.Colormap = cols_14mods;
h1.Position(1) = .54;                  % manual tweak from original code
h1.Position(3) = .02;

% ---- subplot (2): main matrix -------------------------------------------
h2 = subplot(2,2,2);
imagesc(M(ord,ord), clim);
hold on;
box_lines(h2, ci_sorted, n);           % draw thin module borders
h2.Colormap = bwr;
set(h2,'XTick',[],'YTick',[],'TickLength',[0 0]);

% ---- subplot (3): empty placeholder (kept from original) -----------------
subplot(2,2,3,'Position',[0 0 0 0]);   % invisible – for layout consistency

% ---- subplot (4): module colour bar at the bottom ------------------------
h4 = subplot(2,2,4);
imagesc(sort(ci)');                    % a 1‑D colour stripe (horizontal)
set(h4,'YTick',[], ...
        'XTick',module_centres(ci_sorted), ...
        'XTickLabel',mod_names, ...
        'XTickLabelRotation',90, ...
        'FontSize',18,'TickLength',[0 0]);
h4.Colormap = cols_14mods;
h4.Position(2) = .555;
h4.Position(4) = .02;

% -------------------------------------------------------------------------
set(gcf,'Renderer','painters');
end
% ====================== helper functions =================================
function centres = module_centres(ci_sorted)
% Return centre positions of each contiguous block in ci_sorted
bound = [0 ; find(diff(ci_sorted)) ; numel(ci_sorted)];
centres = (bound(1:end-1) + bound(2:end))/2;
end
% -------------------------------------------------------------------------
function box_lines(ax, ci_sorted, n)
% Draw thin black lines at module boundaries
edges = find(diff(ci_sorted));
for e = edges'
    line(ax,[0.5 n+0.5],[e+0.5 e+0.5],'Color','k','LineWidth',0.5);
    line(ax,[e+0.5 e+0.5],[0.5 n+0.5],'Color','k','LineWidth',0.5);
end
end

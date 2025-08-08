function [all_betas all_tstats all_ps]=coupling_parameters(grad_slopes,grad_slope_deltas,grad_slope_2deltas,demean)
% coupling_parameters.m
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

%   [all_betas]=coupling_parameters(grad_slopes,grad_slope_deltas,grad_slope_2deltas)
%   returns the coupling parameters based on a linear model of a given
%   dimension's second derivative (grad_slope_2deltas) as a function of all
%   gradients' timeseries (grad_slopes) and first derivative (grad_slope_deltas)

n_comps=size(grad_slopes,2);
all_betas=[];
all_tstats=[];
all_ps=[];
for i=1:n_comps
    cur_comp=i;
    cur_other_comps=setdiff(1:n_comps,cur_comp);
    cur_X=[];
    cur_y=grad_slope_2deltas(:,cur_comp);
    for j=1:n_comps
        cur_comp_2=j;
        cur_X=[cur_X grad_slopes(:,cur_comp_2) grad_slope_deltas(:,cur_comp_2)];
    end
    if demean
        cur_y=cur_y-mean(cur_y);
        cur_X=cur_X-mean(cur_X);
    end
    cur_mdl=fitlm(cur_X,cur_y);
    %B=robustfit(cur_X,cur_y);
    cur_betas=cur_mdl.Coefficients.Estimate;
    %cur_betas=B;
    cur_tstats=cur_mdl.Coefficients.tStat;
    cur_ps=cur_mdl.Coefficients.pValue;
    all_betas=[all_betas;cur_betas'];
    all_tstats=[all_tstats;cur_tstats'];
    all_ps=[all_ps;cur_ps'];
end
end
function [grad_slope_deltas,grad_slope_2deltas]=latent_derivatives(components_pca,n_subjs,n_vols_per_scan)
% latent_derivatives.m
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

%   [grad_slope_deltas,grad_slope_2deltas]=latent_derivatives(components_pca,n_subjs,n_vols_per_scan) 
%   returns the gradient timeseries first derivative and second derivative
%   based on the gradient timeseries (components_pca), the number of
%   subjects (n_subjs), and the number of fmri volumes per scan (n_vols_per_scan).

n_vols=size(components_pca,1);
n_comps=size(components_pca,2);
grad_slopes=[];
grad_slope_deltas=[];
grad_slope_2deltas=[];
for i=1:n_subjs % calculate derivatives separately for each subject
    cur_offset=((i-1)*n_vols_per_scan)+1;
    cur_inds=cur_offset:cur_offset+(n_vols_per_scan-1);
    grad_slopes=[grad_slopes;components_pca(cur_inds,1:n_comps)];
    temp_grads=zeros(length(cur_inds),n_comps);
    temp_2grads=zeros(length(cur_inds),n_comps);
    for j=1:n_comps % calculate derivatives separately for dimension
        temp_grads(:,j)=gradient(components_pca(cur_inds,j));
        temp_2grads(:,j)=gradient(gradient(components_pca(cur_inds,j)));
    end
    grad_slope_deltas=[grad_slope_deltas;temp_grads];
    grad_slope_2deltas=[grad_slope_2deltas;temp_2grads];
end

end
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

% ========================================================================
% REQUIRED EXTERNAL VARIABLES
% ------------------------------------------------------------------------
% ts_all                 : concatenated BOLD time series [sum(ts_lengths) × n_regions]
%                          Rows stacked by diagnosis in order: {ad, bv, cbs, nfv, sv, hc}.
% ts_lengths             : row counts per group used to split ts_all [6 × 1], order {ad, bv, cbs, nfv, sv, hc}.
% subj_tp_vols_keep2     : row indices (or logical mask) into ts_all selecting the independent-control
%                          cohort used to fit the PCA basis ("group 2"). [vector]
% keep_str_fc            : indices (or logical mask) of scans to include when vectorizing gradient
%                          covariance/FC (maps to 3rd dim of grad_covs_all). [vector of scan indices]

% scanners_ad            : scanner code per AD scan (1 = Trio, 2 = Prisma). [n_scans_ad × 1]
% scanners_bv            : scanner code per bvFTD scan (1 = Trio, 2 = Prisma). [n_scans_bv × 1]
% scanners_cbs           : scanner code per CBS scan (1 = Trio, 2 = Prisma). [n_scans_cbs × 1]
% scanners_nfv           : scanner code per nfvPPA scan (1 = Trio, 2 = Prisma). [n_scans_nfv × 1]
% scanners_sv            : scanner code per svPPA scan (1 = Trio, 2 = Prisma). [n_scans_sv × 1]
% scanners_hc            : scanner code per HC scan (1 = Trio, 2 = Prisma). [n_scans_hc × 1]

% dxs_all                : diagnostic label per scan, coded 1..6 in order {ad, bv, cbs, nfv, sv, hc}. [n_subjs × 1]
% age_all                : age at scan (years) aligned to dxs_all/keep_str_fc. [n_subjs × 1]
% sex_all_tps            : sex code aligned to dxs_all/keep_str_fc (e.g., 1/2). [n_subjs × 1]

% atrophy                : unharmonized atrophy W-scores [n_subjs × n_regions] (ROI order matches ts_all).
% scanners_all           : scanner code per scan for ComBat batch (1 = Trio, 2 = Prisma). [n_subjs × 1]
%
% ========================================================================
% REQUIRED EXTERNAL (USER‑SUPPLIED) FUNCTIONS
% ------------------------------------------------------------------------
% combat                – ComBat batch normalization (see https://github.com/Jfortin1/ComBatHarmonization/tree/master/Matlab)

restoredefaultpath; rehash toolboxcache; cd(fileparts(mfilename('fullpath'))); addpath(genpath(pwd));

% get gradients for independent controls
load('adftd_structure_function_data.mat');
parts = mat2cell(ts_all, ts_lengths(:), size(ts_all,2));  % split by rows
[ts_ad, ts_bv, ts_cbs, ts_nfv, ts_sv, ts_hc] = parts{:};  % unpack

keep_regions=1:246;
[roi_comp_slopes_fc_in_grp2 codes_fc_in_grp2 x x explained_in_grp2]=pca(ts_all(subj_tp_vols_keep2,keep_regions));

% project patients into control space
% see https://www.mathworks.com/matlabcentral/answers/53259-how-to-project-a-new-point-to-pca-new-basis
codes_fc_in_all_proj_grp2=(ts_all(:,keep_regions)-mean(ts_all(subj_tp_vols_keep2,keep_regions)))/roi_comp_slopes_fc_in_grp2';

clear codes_all codes_fc_in codes_fc_in_grp1 codes_fc_in_grp2 
clear cur_dx_ts cur_dx_codes cur_dx_codes_proj_grp2 fcmats_all_flat
clear ts_all

% get gradient covariance
n_comps=246;
clear grad_covs_all coupling_all grad_covs_all_proj_grp2
count=1;
dxs={'ad' 'bv' 'cbs' 'nfv' 'sv' 'hc'};
ts_lengths=[length(ts_ad) length(ts_bv) length(ts_cbs) length(ts_nfv) length(ts_sv) length(ts_hc)];
get_coupling=true;
for i=1:6
    cur_dx_scanners=eval(sprintf('scanners_%s',dxs{i}));
    cur_dx_n_scans=length(cur_dx_scanners);
    cur_dx_ts=eval(sprintf('ts_%s',dxs{i}));
    cur_dx_codes_proj_grp2=codes_fc_in_all_proj_grp2((sum(ts_lengths(1:i-1)))+1:sum(ts_lengths(1:i)),:);
    start_n=1;
    for j=1:cur_dx_n_scans
        if cur_dx_scanners(j)==2 % prisma
            n_vols=555;
        else % trio
            n_vols=235;
        end
        stop_n=start_n+n_vols-1;

        cur_codes=cur_dx_codes_proj_grp2(start_n:stop_n,:);

        grad_covs_all(:,:,count)=cov(cur_codes(:,1:n_comps));
        cur_grad_slopes=cur_codes(:,1:6);

        subj_n_vols(count)=n_vols;
        disp(count)
        count=count+1;
        start_n=start_n+n_vols;
    end
end

% vectorize gradient covariance matrices
count=1;
n_comps=246;
keep_inds=keep_str_fc;
clear pls_grad_cov pls_grad_pairs pls_grad_cov_proj_grp2
for i=1:n_comps
    for j=i:n_comps
        pls_grad_cov(:,count)=squeeze(grad_covs_all(i,j,keep_inds));
        pls_grad_cov(:,count)=squeeze(grad_covs_all(i,j,keep_inds));
        pls_grad_pairs(count,1)=i;
        pls_grad_pairs(count,2)=j;
        count=count+1;
    end
end

% gradient region weights
grad_weights=roi_comp_slopes_fc_in_grp2;

% combat to remove batch effects of scanner for atrophy and gradient covariance/fc
n_grads=6;
n_subjs=length(keep_inds);
n_pts=length(find(dxs_all(keep_inds)<6));
n_regions=length(keep_regions);
combat_covs=dummyvar(dxs_all(keep_inds));
combat_covs=combat_covs(:,2:end);
grp_cov=zeros(n_subjs,1);grp_cov(find(dxs_all(keep_inds)<6))=1;
combat_covs=[grp_cov age_all(keep_inds) sex_all_tps(keep_inds)-1];
atrophy_300_harmonized=combat(atrophy(keep_inds,1:n_regions)',scanners_all(keep_inds),combat_covs,1)';
pls_grad_cov_300_harmonized=combat(pls_grad_cov',scanners_all(keep_inds),combat_covs,1)';
clear grad_covs_300_harmonized pls_grad_cov_n pls_grad_pairs_n
count=1;
for i=1:n_grads
    for j=i:n_grads
        cur_ind=intersect(find(pls_grad_pairs(:,1)==i),find(pls_grad_pairs(:,2)==j));
        grad_covs_300_harmonized(i,j,:)=pls_grad_cov_300_harmonized(:,cur_ind);
        grad_covs_300_harmonized(j,i,:)=pls_grad_cov_300_harmonized(:,cur_ind);
        pls_grad_cov_n(:,count)=squeeze(grad_covs_300_harmonized(i,j,:))';
        pls_grad_pairs_n(count,1)=i;
        pls_grad_pairs_n(count,2)=j;
        count=count+1;
    end
end

% reconstruct FC matrix from gradient covariance after combat
k1=sub2ind([n_comps n_comps],pls_grad_pairs(:,1),pls_grad_pairs(:,2));
k2=sub2ind([n_comps n_comps],pls_grad_pairs(:,2),pls_grad_pairs(:,1));
fc_mats_via_codes=zeros(n_regions,n_regions,n_subjs);
grad_covs_all_harmonized=zeros(n_comps,n_comps,n_subjs);
for i=1:n_subjs
    cur_grad_covs=zeros(n_comps,n_comps);
    cur_grad_covs(k1)=pls_grad_cov_300_harmonized(i,:);
    cur_grad_covs(k2)=pls_grad_cov_300_harmonized(i,:);
    grad_covs_all_harmonized(:,:,i)=cur_grad_covs;
    var1_all=diag(grad_weights(1:n_regions,1:n_comps)*cur_grad_covs*grad_weights(1:n_regions,1:n_comps)');
    cov1_all=grad_weights(1:n_regions,1:n_comps)*cur_grad_covs*grad_weights(1:n_regions,1:n_comps)';
    var1_2_all=sqrt(var1_all*var1_all');
    fc_mats_via_codes(:,:,i)=cov1_all./var1_2_all;
    disp(i)
end

% vectorize FC for PLS
count=1;
clear pls_fc pls_fc_pairs pls_cov
pls_fc=[];
pls_fc_pairs=[];
for i=1:length(keep_regions)-1
    for j=i+1:length(keep_regions)
        pls_fc(:,count)=squeeze(fc_mats_via_codes(i,j,:));
        pls_fc_pairs(count,1)=i;
        pls_fc_pairs(count,2)=j;
        count=count+1;
    end
end
fcmats_mean=mean(pls_fc,2);
% linear indices for upper triangle matrix edges
keep_pls_fc=sub2ind([n_regions n_regions],pls_fc_pairs(:,1),pls_fc_pairs(:,2));
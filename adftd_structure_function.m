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
% atrophy_300_harmonized : atrophy W‑scores (ComBat‑harmonized) [scans × 246 ROIs]
% pls_fc                : vectorized FC upper‑triangle  (ComBat‑harmonized) [scans × 30135 edges]
% keep_pls_fc           : linear indices of the upper‑triangle edges (maps Yloadings → 246×246)
% fc_mats_via_codes     : full 246×246 FC matrices for each scan
% keep_str_fc           : indices of scans with usable structural‑FC data
% dxs_all               : integer diagnosis label (1–6) for each scan
% keep_ad_sub           : logical / index vector for Alzheimer’s disease scans (1)
% keep_bv_sub           : “             ” behavioural‑variant FTD (2)
% keep_cbs_sub          : “             ” corticobasal syndrome (3)
% keep_nfv_sub          : "             ” nonfluent variant PPA (4)
% keep_sv_sub           : “             ” semantic variant PPA (5)
% keep_hc_sub           : “             ” healthy controls (6)
% ci_15mods_v2          : community / network ID (1–15) per ROI (Brown et al. 2019)
% grad_covs_all_harmonized : 6×6×N tensor of gradient covariances (post‑ComBat)
% pls_grad_cov_n        : vectorized gradient covariances [scans × 21 pairs]
% grad_weights (roi_comp_slopes_fc_in_grp2)
%                       : ROI weights mapping gradient scores → FC (246 × nComp)
% scanners_all          : scanner ID per scan (1 = Trio, 2 = Prisma)
% codes_fc_in_all_proj_grp2 : concatenated gradient‑slope time‑series [timepts × 6]
% subj_tp_vols          : mapping from each row of codes_fc_in_all_proj_grp2 → scan index
% pls_fc_pairs          : 30135 × 2 list of ROI pairs matching keep_pls_fc order
% pls_grad_pairs_n      : 21 × 2 list of gradient‑pair indices matching pls_grad_cov_n
%
% ========================================================================
% REQUIRED EXTERNAL (USER‑SUPPLIED) FUNCTIONS
% ------------------------------------------------------------------------
% plot_matrix                – helper to visualise FC / loading matrices
% plotcorr                   – scatter‑plot with correlation line
% latent_derivatives         – compute first/second derivatives of gradient slopes
% coupling_parameters        – estimate coupled‑oscillator parameters
% coupled_oscillator_function– simulate gradient dynamics from coupling matrix
% circ_mean                  – circular mean (e.g., CircStat toolbox)
% wrapTo2Pi                  – wrap angles to [0, 2π] (Mapping Toolbox; built‑in)
% ========================================================================

%load('adftd_structure_function_data.mat');
n_regions=246; % brainnetome cortical and subcortical regions
keep_regions=1:246;
n_grads=6; % keep first six gradients
% first 321 scans are baseline scans for subjects 1-321; scans 321-446 are
% follow-up scans for subset of subjects with longitudinal scans
baseline_inds=1:321;
% gradient covariance matrices
grad_covs_all_mean=mean(grad_covs_all_harmonized(:,:,baseline_inds),3);
% region gradient weights
grad_weights=roi_comp_slopes_fc_in_grp2;

%% PLS on atrophy [321 x 246] versus functional connectivity [321 x 30135] (Figure 1)
first_321=baseline_inds;
[Xloadings,Yloadings,Xscores,Yscores,betapls,plspctvar,mse,plsstats]=plsregress(atrophy_300_harmonized(first_321,keep_regions),pls_fc(first_321,:),10);

% flip signs and weights for component 1 so high scores = high atrophy
Xscores(:,1)=-Xscores(:,1);
Yscores(:,1)=-Yscores(:,1);
Xloadings(:,1)=-Xloadings(:,1);
Yloadings(:,1)=-Yloadings(:,1);
plsstats.W(:,1)=-plsstats.W(:,1);

% plot PLS component score correlations
figure;
for comp_num=1:5
    subplot(2,3,comp_num)
    plotcorr(Xscores(:,comp_num),Yscores(:,comp_num))
    title(sprintf('PLS comp %d',comp_num))
end

figure;
plot_reg_line=true;
for comp_num=1:5
    cur_data=[Xscores(:,comp_num) Yscores(:,comp_num)];
    subplot(2,3,comp_num)
    scatter(cur_data(:,1),cur_data(:,2),50)
    hold on
    scatter(cur_data(keep_ad_sub,1),cur_data(keep_ad_sub,2),50,'filled','red')
    scatter(cur_data(keep_bv_sub,1),cur_data(keep_bv_sub,2),50,'filled','blue')
    scatter(cur_data(keep_cbs_sub,1),cur_data(keep_cbs_sub,2),50,'filled','yellow')
    scatter(cur_data(keep_nfv_sub,1),cur_data(keep_nfv_sub,2),50,'filled','magenta')
    scatter(cur_data(keep_sv_sub,1),cur_data(keep_sv_sub,2),50,'filled','cyan')
    scatter(cur_data(keep_hc_sub,1),cur_data(keep_hc_sub,2),50,'filled','green')
    title(sprintf('PLS comp %d',comp_num))
    xlabel('atrophy score')
    ylabel('FC score')
    if plot_reg_line
        y=Yscores(:,comp_num);
        x=Xscores(:,comp_num);
        [b,bint,r,rint,stats]=regress(y,[ones(length(x),1) x]);
        % plot regression line
        x_span=max(x)-min(x);
        x_buffer=x_span*.1;
        xreg=linspace(min(x)-x_buffer,max(x)+x_buffer,10);
        yreg=b(1)+(b(2).*xreg);
        regline=plot(xreg , yreg,'-');
        set(regline,'LineWidth',2,'Color',[0 0 0]);
        axis tight
    end
end

% plot PLS component FC loading matrices
clear pls_yloadings_mats_all
for comp_num=1:5
    pls_yloadings_mat=zeros(n_regions,n_regions);
    for i=1:((n_regions*(n_regions-1))/2)
        pls_yloadings_mat(pls_fc_pairs(i,1),pls_fc_pairs(i,2))=Yloadings(i,comp_num);
    end
    pls_yloadings_mat=pls_yloadings_mat+pls_yloadings_mat';
    if comp_num==1
        pls_yloadings_mat=-pls_yloadings_mat;
    end
    pls_yloadings_mats_all(:,:,comp_num)=pls_yloadings_mat;
end

for comp_num=1:5
    if comp_num==1
        fcmat=-pls_yloadings_mats_all(:,:,comp_num);
    else
        fcmat=pls_yloadings_mats_all(:,:,comp_num);
    end
    plot_matrix(fcmat,ci_15mods_v2,[-1 1])
end

%% Multidimensional scaling of structural components (Figure 2)
options = statset('MaxIter',5000);
D=pdist(Xscores(:,1:3),"seuclidean");
reduction = mdscale(D,2,'Options',options);
cur_data=reduction;
figure('position',[0 0 800 800])
dot_size=400;
scatter(cur_data(:,1),cur_data(:,2),dot_size)
hold on
scatter(cur_data(keep_hc_sub,1),cur_data(keep_hc_sub,2),dot_size,'filled','green','MarkerEdgeColor',[0 0 0])
scatter(cur_data(keep_ad_sub,1),cur_data(keep_ad_sub,2),dot_size,'filled','red','MarkerEdgeColor',[0 0 0])
scatter(cur_data(keep_bv_sub,1),cur_data(keep_bv_sub,2),dot_size,'filled','blue','MarkerEdgeColor',[0 0 0])
scatter(cur_data(keep_cbs_sub,1),cur_data(keep_cbs_sub,2),dot_size,'filled','yellow','MarkerEdgeColor',[0 0 0])
scatter(cur_data(keep_nfv_sub,1),cur_data(keep_nfv_sub,2),dot_size,'filled','magenta','MarkerEdgeColor',[0 0 0])
scatter(cur_data(keep_sv_sub,1),cur_data(keep_sv_sub,2),dot_size,'filled','cyan','MarkerEdgeColor',[0 0 0])

dx_colors={'red','blue','yellow','magenta','cyan','green'};
if true
clear dx_centroids dx_centroids_3d
for i=1:6
    cur_inds=find(dxs_all(keep_str_fc)==i);
    cur_dx_centroid=mean(cur_data(cur_inds,:));
    cur_dx_centroid_3d=mean(Xscores(cur_inds,1:3));
    dx_centroids(i,:)=cur_dx_centroid;
    dx_centroids_3d(i,:)=cur_dx_centroid_3d;
    scatter(cur_dx_centroid(1),cur_dx_centroid(2),1000,'filled',dx_colors{i},'MarkerEdgeColor','black','LineWidth',6)
end
end

set(gca,'FontSize',24);
xlabel('Dim 1')
ylabel('Dim 2')
f = gcf;
f.Renderer = 'painters';
Tight = get(gca, 'TightInset');
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
set(gca, 'Position', NewPos);


%% Syndrome-specific FC matrices (Figure 2)
% linear discriminant analysis to get patients with typical atrophy pattern for each
% syndrome
lda_data=[Xscores(:,1:3)];
lda_grp=dxs_all(keep_str_fc);
rng('default') % For reproducibility
[class,err,posterior,logp,coeff]=classify(lda_data,lda_data,lda_grp);
figure;
cm = confusionchart(lda_grp,class);
cur_inds_correct_all=[];
for i=1:6
    cur_inds=find(dxs_all(keep_str_fc)==i);
    cur_inds_correct=cur_inds(find((class(cur_inds)-lda_grp(cur_inds))==0));
    cur_inds_correct_all=[cur_inds_correct_all;cur_inds_correct];
end
cur_inds_incorrect_all=setdiff(1:321,cur_inds_correct_all);
classify_correct=zeros(321,1);classify_correct(cur_inds_correct_all)=1;

[B,dev,stats] = mnrfit(lda_data,lda_grp);

% fc matrix differences for each syndrome vs cognitively normal subjects
clear mean_fc_k1_sub_k2
n_clusters=6;
for i=1:n_clusters
    ref_cluster=6;
    cur_k1=intersect(find(dxs_all(keep_str_fc)==i),cur_inds_correct_all);
    cur_k2=intersect(find(dxs_all(keep_str_fc)==ref_cluster),cur_inds_correct_all);
    mean_fc_k1=mean(fc_mats_via_codes(:,:,cur_k1),3);
    mean_fc_k2=mean(fc_mats_via_codes(:,:,cur_k2),3);
    mean_fc_k1_all(:,:,i)=mean_fc_k1;
    mean_fc_k1_sub_k2(:,:,i)=mean_fc_k1-mean_fc_k2;
    if false
        for j=1:n_regions
            for k=1:n_regions
                [h p ci stats]=ttest2(squeeze(fc_mats_via_codes(j,k,cur_k1)),squeeze(fc_mats_via_codes(j,k,cur_k2)));
                t_fc_k1_sub_k2(j,k,i)=stats.tstat;
            end
        end
    end
    disp(i)
end

% reconstruct syndrome fc patterns from function components 1-3 (Yscores)
clear mean_fc_recon_clusters
n_clusters=6;
for i=1:n_clusters
    cur_k=intersect(find(dxs_all(keep_str_fc)==i),cur_inds_correct_all);
    cur_fc=Yscores(cur_k,1:3)*Yloadings(:,1:3)';
    cur_fc_mean=mean(cur_fc);
    mean_fc_recon=zeros(n_regions);
    mean_fc_recon(keep_pls_fc)=cur_fc_mean;
    mean_fc_recon=mean_fc_recon+mean_fc_recon';
    mean_fc_recon_clusters(:,:,i)=mean_fc_recon;
end

% fc matrix reconstructed differences for each syndrome vs cognitively normal subjects
clear mean_fc_recon_k1_sub_k2
for i=1:n_clusters
    ref_cluster=6;
    mean_fc_recon_k1_sub_k2(:,:,i)=mean_fc_recon_clusters(:,:,i)-mean_fc_recon_clusters(:,:,ref_cluster);
end

% combined matrix with actual fc on upper triangle, reconstructed fc on
k_ut=find(triu(ones(246),1));
k_lt=find(tril(ones(246),-1));

syndrome_abbrevs={'ad','bv','cbs','nfv','sv'};
[sci ici]=sort(ci_15mods_v2);
for cur_ind=1:5
    cur_mat1=mean_fc_k1_sub_k2(:,:,cur_ind);

    % reconstructed fc
    cur_mat2=mean_fc_recon_k1_sub_k2(:,:,cur_ind);
    cur_mat2=cur_mat2./max(abs(cur_mat2(:)));
    cur_mat2=cur_mat2.*max(abs(cur_mat1(:)));

    cur_mean_fc=zeros(n_regions);
    for j=1:(n_regions-1)
        cur_ind1=ici(j);
        for k=(j+1):n_regions
            cur_ind2=ici(k);
            [s2 i2]=sort([cur_ind1 cur_ind2]);
            cur_mean_fc(j,k)=cur_mat1(cur_ind1,cur_ind2);
            cur_mean_fc(k,j)=cur_mat2(cur_ind1,cur_ind2);
        end
    end

    mean_fc_dx_all(:,:,cur_ind)=cur_mean_fc;

    fcmat=cur_mean_fc;
    % note: matrices are resorted by modular order, so pass sorted module
    % sequence, rather than sequence based on original region ordering
    plot_matrix(fcmat,sci,[-.2 .2]);
end


%% Correlation between PLS functional components and gradient covariance (Figure 4)
% variance in Yscores explained by gradients
pls_grad_cov_n_z=zscore(pls_grad_cov_n(baseline_inds,:));
z_yscores=zscore(Yscores(baseline_inds,1:3));
for comp_num=1:3
    stats=regstats(z_yscores(baseline_inds,comp_num),pls_grad_cov_n_z);
    n_tests=(n_grads*(n_grads+1))/2;
    k=find(stats.tstat.pval(2:end)<(.05/n_tests));
    ts=stats.tstat.t(2:end);ps=stats.tstat.pval(2:end);
    [s1 i1]=sort(abs(ts),'descend');
    
    % First line: r‑squared and count of significant gradient pairs
    fprintf('r-squared: %.4f, # of significant gradient pairs: %d\n', ...
            stats.rsquare, length(k));
    
    % Second line: column headings
    fprintf('grad A\tgrad B\tt-stat\tp-value\n');
    
    % Then print each row of [grad A, grad B, t, p]
    for idx = 1:length(k)
        i = i1(idx);
        gA   = pls_grad_pairs_n(i,1);
        gB   = pls_grad_pairs_n(i,2);
        tval = ts(i);
        pval = ps(i);
        fprintf('%4d\t%4d\t%8.3f\t%8.4f\n', gA, gB, tval, pval);
    end
end

% variance percent per predictor
% helpful: https://www.mathworks.com/help/stats/coefficient-of-determination-r-squared.html
% this approach for estimating per-term variance is based on:
% https://www.researchgate.net/publication/306347340_A_Natural_Decomposition_of_R2_in_Multiple_Linear_Regression
for comp_num=1:3
    mdl=fitlm(pls_grad_cov_n_z,z_yscores(:,comp_num));
    clear term_rsqs
    for i=1:21
        c=cov(pls_grad_cov_n_z(:,i),mdl.Fitted);
        term_rsqs(i)=c(1,2)*mdl.Coefficients.Estimate(i+1);
    end
    figure('position',[0 0 800 800])
    [srsq irsq]=sort(term_rsqs(1:21),'descend');
    bar([srsq'],'FaceColor','k','EdgeColor','k','LineWidth',0.01)
    max(cumsum(srsq))
    ticks_labels={};
    for i=1:length(irsq)
        cur_grads=pls_grad_pairs_n(irsq(i),:);
        ticks_labels{i}=sprintf('%d-%d',cur_grads(1),cur_grads(2));
    end
    xticks(1:21)
    xticklabels(ticks_labels)
    xtickangle(90)
    ylim([0 max(srsq)])
    set(gca,'FontSize',36);
    Tight = get(gca, 'TightInset');
    NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-Tight(4)]; %New plot position [X Y W H]
    set(gca, 'Position', NewPos);
end

% show subject A - subject B partial difference fc matrices based on one
% specific gradient pair
n_comps=246;
grad_covs_all_mean=mean(grad_covs_all_harmonized(:,:,baseline_inds),3);

subj_mat=[
    272 202;
    272 202;
    2 199;
    124 297];
grad_mat=[
    1 1;
    1 4;
    1 2;
    2 2];
display_thresholds=[
    -.3 .3;
    -.25 .25;
    -.3 .3;
    -.15 .15];

% plot four examples from paper
for i=1:4
    cur_subj_1=subj_mat(i,1);
    cur_subj_2=subj_mat(i,2);
    grad1=grad_mat(i,1);
    grad2=grad_mat(i,2);

    % get subject A fc matrix
    % 1) start with group mean gradient covariance matrix
    % 2) for specific gradient or pair of gradients, replace the mean
    % variance/covariance values with those from this subject (deviation from
    % the mean)
    % 3) get resultant fc matrix
    grad_covs_all_subj1=grad_covs_all_mean;
    grad_covs_all_subj1(grad1,grad2)=grad_covs_all_harmonized(grad1,grad2,cur_subj_1);
    grad_covs_all_subj1(grad2,grad1)=grad_covs_all_harmonized(grad2,grad1,cur_subj_1);
    
    var1_all=diag(grad_weights(1:n_regions,1:n_comps)*grad_covs_all_subj1*grad_weights(1:n_regions,1:n_comps)');
    cov1_all=grad_weights(1:n_regions,1:n_comps)*grad_covs_all_subj1*grad_weights(1:n_regions,1:n_comps)';
    var1_2_all=sqrt(var1_all*var1_all');
    fc_from_grad_covs_subj1=cov1_all./var1_2_all;
    
    % get subject B fc matrix
    grad_covs_all_subj2=grad_covs_all_mean;
    grad_covs_all_subj2(grad1,grad2)=grad_covs_all_harmonized(grad1,grad2,cur_subj_2);
    grad_covs_all_subj2(grad2,grad1)=grad_covs_all_harmonized(grad2,grad1,cur_subj_2);
    
    var1_all=diag(grad_weights(1:n_regions,1:n_comps)*grad_covs_all_subj2*grad_weights(1:n_regions,1:n_comps)');
    cov1_all=grad_weights(1:n_regions,1:n_comps)*grad_covs_all_subj2*grad_weights(1:n_regions,1:n_comps)';
    var1_2_all=sqrt(var1_all*var1_all');
    fc_from_grad_covs_subj2=cov1_all./var1_2_all;
    
    % compute and plot difference matrix
    %fc_from_grad_covs_dif=fc_from_grad_covs_subj2-fc_from_grad_covs_subj1;
    fc_from_grad_covs_dif=fc_from_grad_covs_subj1-fc_from_grad_covs_subj2;
    fcmat=fc_from_grad_covs_dif;
    plot_matrix(fcmat,ci_15mods_v2,display_thresholds(i,:))
end


%% Coupled oscillator analysis, comparison with gradient covariance (Figure 5/6)
window_size=1;
step_size=1;
keep_grads=1:6;
clear coupling_windowed Yscores_windowed grad_covs_real_windowed grad_covs_pred_windowed fcmat_real_windowed fcmat_pred_windowed v_windowed d_windowed freqs_windowed ks_mean_windowed
count=1;
pred_interval=235;
for i=1:step_size:length(baseline_inds)-(window_size-1)
    cur_inds=baseline_inds(i:i+(window_size-1)); % 321 indexed
    if length(cur_inds)==1
        cur_scanner=scanners_all(keep_str_fc(cur_inds));
        if cur_scanner==1
            pred_interval=235;
        else
            pred_interval=555;
        end
    end

    cur_subj_tp_vols=[];
    for j=1:window_size
        cur_subj_tp_vols=[cur_subj_tp_vols;find(subj_tp_vols==(keep_str_fc(cur_inds(j))))];
    end
    % concatenate gradient timeseries for current group of subjects
    cur_grads=codes_fc_in_all_proj_grp2(cur_subj_tp_vols,keep_grads);
    
    % compute real coupling parameters, gradient covariance, Yscores, and
    % functional connectivity
    n_vols=size(cur_grads,1);
    cur_grad_slopes=cur_grads;
    [cur_grad_slope_deltas cur_grad_slope_2deltas]=latent_derivatives(cur_grad_slopes,1,n_vols);
    [cur_betas]=coupling_parameters(cur_grad_slopes,cur_grad_slope_deltas,cur_grad_slope_2deltas,true);
    grad_covs_real_windowed(:,:,count)=cov(cur_grad_slopes);
    coupling_windowed(:,:,count)=cur_betas;
    
    % eigendecomposition of coupling parameters
    A=zeros(n_grads*2,n_grads*2);
    for j=1:n_grads
        i1=(j*2)-1;
        i2=j*2;
        A(i1,i2)=1;
        A(i2,:)=cur_betas(j,2:((n_grads*2)+1));
    end
    [v d]=eig(A);
    v_windowed(:,:,count)=v;
    d_windowed(:,:,count)=d;
    
    % eigenmode frequency measures
    ang_freqs=diag(imag(d(1:2:end,1:2:end))); % number of 2pi radians covered in one data timestep
    nat_periods=(2*pi)./ang_freqs; % number of timesteps to cover one full cycle
    timestep_length=2; % length in seconds of one data timestep (TR)
    freqs=1./(nat_periods*timestep_length); % number of cycles per second (Hz)
    freqs_windowed(count,:)=freqs';
    
    if length(cur_inds)>1
        Yscores_windowed(count,:)=mean(Yscores(cur_inds,:));
    else
        Yscores_windowed(count,:)=Yscores(cur_inds,:);
    end
    cur_real_region_ts=cur_grad_slopes*grad_weights(1:n_regions,keep_grads)';
    fcmat_real_windowed(:,:,count)=corr(cur_real_region_ts);
    
    % simulate gradient timeseries using current coupling parameters
    clear cur_grad_covs_pred_windowed cur_fcmat_pred_windowed cur_ks

    [pos_vel_mat t ks_keep]=coupled_oscillator_function(cur_betas,cur_grad_slopes,cur_grad_slope_deltas,grad_weights,pred_interval);

    for j=1:30
        start=j;
        ic1=cur_grad_slopes(start,keep_grads);
        ic2=cur_grad_slope_deltas(start,keep_grads);
        ics=[];
        for cur_grad=1:n_grads
            ics=[ics;ic1(cur_grad)];
            ics=[ics;ic2(cur_grad)];
        end
        ks=linsolve(pos_vel_mat(:,ks_keep,1),ics);
        
        cur_ks(:,j)=ks;
        clear pos_vel
        for cur_t=1:length(t)
            pos_vel(:,cur_t)=pos_vel_mat(:,ks_keep,cur_t)*ks;
        end
        cur_y_pos_all=pos_vel(1:2:end,:);
        cur_grad_covs_pred_windowed(:,:,j)=cov(cur_y_pos_all');
        cur_sim_region_ts=(cur_y_pos_all)'*grad_weights(1:n_regions,keep_grads)';
        fcmat=corr(cur_sim_region_ts);
        cur_fcmat_pred_windowed(:,:,j)=fcmat;
    end

    grad_covs_pred_windowed(:,:,count)=mean(cur_grad_covs_pred_windowed,3);
    fcmat_pred_windowed(:,:,count)=mean(cur_fcmat_pred_windowed,3);
    ks_mean_windowed(:,count)=mean(cur_ks,2);
    
    disp(count)
    count=count+1;
end


%% Get pure eigenmode fc matrices (Figure 5)
window_size=5;
step_size=1;
k1=find(scanners_all(keep_str_fc)==1); % trio only
[s1 i1]=sort(Yscores(:,1));
i12=intersect(i1,k1,'stable');
keep_grads=1:6;
n_grads=length(keep_grads);

% start index for window
window_starts=[1 168];
for i=1:2
    cur_window_start=window_starts(i);
    low_atrophy=false;
    if cur_window_start==1
        low_atrophy=true;
    end

    cur_inds=i12(cur_window_start:cur_window_start+(window_size-1)); % 321 indexed
    cur_subj_tp_vols=[];
    for j=1:window_size
        cur_subj_tp_vols=[cur_subj_tp_vols;find(subj_tp_vols==(keep_str_fc(cur_inds(j))))];
    end
    % concatenate gradient timeseries for current group of subjects
    cur_grad_slopes=codes_fc_in_all_proj_grp2(cur_subj_tp_vols,keep_grads);
    % compute coupling parameters
    n_vols=size(cur_grad_slopes,1);
    [cur_grad_slope_deltas cur_grad_slope_2deltas]=latent_derivatives(cur_grad_slopes,1,n_vols);
    [cur_betas]=coupling_parameters(cur_grad_slopes,cur_grad_slope_deltas,cur_grad_slope_2deltas,true);

    cur_pred_interval=235;
    start_condition=1;
    [pos_vel_mat t ks_keep]=coupled_oscillator_function(cur_betas,cur_grad_slopes,cur_grad_slope_deltas,grad_weights,cur_pred_interval,start_condition);

    clear mode_ts_all
    for mode_num=1:6
        mode_real=(4*(mode_num-1))+1;
        mode_imag=(4*(mode_num-1))+2;
        mode_ts=squeeze(pos_vel_mat(1:2:12,mode_real,:)+pos_vel_mat(1:2:12,mode_imag,:))';
        mode_ts_all(:,:,mode_num)=mode_ts;
    end

    mode_ts_sum=sum(mode_ts_all,3);
    roi_ts_sum=mode_ts_sum*grad_weights(:,1:6)';
    fcmat=corr(roi_ts_sum);

    if low_atrophy
        coupling_windowed_low_atrophy=cur_betas;
        mode_ts_all_low_atrophy=mode_ts_sum;
        pos_vel_mat_low_atrophy=pos_vel_mat;
        fcmat_sim_low_atrophy=fcmat;
    else
        coupling_windowed_high_atrophy=cur_betas;
        mode_ts_all_high_atrophy=mode_ts_sum;
        pos_vel_mat_high_atrophy=pos_vel_mat;
        fcmat_sim_high_atrophy=fcmat;
    end
end

plot_matrix(fcmat_sim_low_atrophy,ci_15mods_v2,[-1 1])
plot_matrix(fcmat_sim_high_atrophy,ci_15mods_v2,[-1 1])
plot_matrix(fcmat_sim_high_atrophy-fcmat_sim_low_atrophy,ci_15mods_v2,[-.5 .5])


% gradient amplitude + gradient pair angles for each eigenmode
clear grad_amps grad_pair_angs grad_pair_angs_mean grad_mags_prods grad_pair_angs_weightmean
n_windows=321;
[all_j all_k]=find(triu(ones(n_grads),1));
grad_amps=zeros(n_windows,n_grads);
for i=1:n_windows
    cur_eigvecs=squeeze(v_windowed(:,:,i));
    cur_eigvals=squeeze(d_windowed(:,:,i));

    % k constants
    % option 1: all ones
    ks=ones(n_grads*2,1);

    % get amplitude for each gradient on each mode
    % eigenvector magnitude represents amplitude of each gradient's fluctuation
    % on a given eigenmode
    clear cur_grad_amps
    amps_real=real(cur_eigvecs(:,1:2:end));
    amps_imag=imag(cur_eigvecs(:,1:2:end));
    for j=1:n_grads*2
        for k=1:n_grads
            cur_real=amps_real(j,k);
            cur_imag=amps_imag(j,k);
            cur_k_real=ks((k*2)-1);
            cur_k_imag=ks(k*2);
            cur_grad_amps(j,k)=sqrt(((cur_k_real*(-cur_real+cur_imag))^2) + ((cur_k_imag*(cur_real+cur_imag))^2));
        end
    end
    cur_grad_amps=cur_grad_amps(1:2:end,:); % trim to only position
    grad_amps(i,:)=sum(cur_grad_amps,2)'; % sum amplitude across modes

    % get angle for each gradient pair on each mode
    cur_grad_angs=wrapTo2Pi(angle(cur_eigvecs(1:2:end,1:2:end)));
    for n=1:length(all_j)
        j=all_j(n);
        k=all_k(n);
        % angle between gradient pairs on each mode
        normdeg=mod([cur_grad_angs(j,:)-cur_grad_angs(k,:)]',2*pi);
        absdiffdeg=min((2*pi)-normdeg,normdeg);
        grad_pair_angs(i,n,:)=absdiffdeg';
        % unweighted average
        grad_pair_angs_mean(i,n)=circ_mean(absdiffdeg);
        
        % weighted average based on gradient magnitude products
        cur_grad_amps_prods=cur_grad_amps(j,:).*cur_grad_amps(k,:);
        grad_pair_angs_weightmean(i,n)=circ_mean(absdiffdeg,cur_grad_amps_prods');
    end
end


%% Correlate gradient covariance with eigenmode magnitude/angle (Supplementary Figure 10)
[all_j all_k]=find(triu(ones(n_grads),0));
[all_j1 all_k1]=find(triu(ones(n_grads),1));
ticks_labels={};
clear grad_cov_sub grad_mag_ang_sub
for i=1:length(all_j)
    cur_grad1=all_j(i);
    cur_grad2=all_k(i);
    
    % either use 'raw' gradient covariance, or post-combat gradient covariance
    % reminder: eigenmodes are computed from 'raw' gradient timeseries before combat
    if false
        cur_grad_cov=squeeze(grad_covs_real_windowed(cur_grad1,cur_grad2,:));
    else
        cur_ind=intersect(find(pls_grad_pairs_n(:,1)==cur_grad1),find(pls_grad_pairs_n(:,2)==cur_grad2));
        cur_grad_cov=pls_grad_cov_n(baseline_inds,cur_ind);
    end

    grad_cov_sub(:,i)=cur_grad_cov;
    if cur_grad1==cur_grad2
        grad_mag_ang_sub(:,i)=grad_amps(:,cur_grad1);
    else
        ind1=find(all_j1==all_j(i));
        ind2=find(all_k1==all_k(i));
        ind3=intersect(ind1,ind2);
        grad_mag_ang_sub(:,i)=grad_pair_angs_weightmean(:,ind3);
    end
    ticks_labels{i}=sprintf('%d-%d',cur_grad1,cur_grad2);
end
[grad_var_mode_corrs grad_var_mode_corr_ps]=corr(grad_cov_sub,grad_mag_ang_sub);


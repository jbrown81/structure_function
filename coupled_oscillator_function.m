function [pos_vel_mat t ks_keep fcmat] = coupled_oscillator_function(coupling_params,grad_slopes,grad_slope_deltas,grad_weights,pred_interval,start_condition)

% coupled_oscillator_model.m
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

% sets up and runs coupled oscillator model based on activity gradients and
% their coupling parameters
% for background, see https://doi.org/10.1016/j.neuroimage.2022.119526

% requirements
% data:
% coupling_params: coupling parameters
% grad_slopes: gradient slope timeseries
% grad_slope_deltas: gradient slope first derivative timeseries
% grad_weights: region weights on each gradient
% functions:
% gradient_ode.m

if nargin < 6
    start_condition=1;
end

% specify parameters
keep_grads=1:6;
%pred_interval=235;
step_size=1;
n_grads=length(keep_grads);
%grad_weights=roi_comp_slopes_fc_in_grp2;

%coupling_params=coupling_windowed(1:n_grads,1:((n_grads*2)+1),1);
%coupling_params=coupling_windowed(1:n_grads,1:((n_grads*2)+1),end);

% set up coupling parameters as system of first order equations in matrix form
% get eigenvalues + eigenvectors
A=zeros(n_grads*2,n_grads*2);
for i=1:n_grads
    i1=(i*2)-1;
    i2=i*2;
    A(i1,i2)=1;
    A(i2,:)=coupling_params(i,2:((n_grads*2)+1));
end
[v d]=eig(A);

% initial conditions vector, two options:
% 1) select a specific timepoint from the actual gradient timeseries data
% 2) "excite" a specific eigenmode
% note that sequence of eigenmodes determined by 'eig' are in random order,
% so resort them from low to high frequency if needed

if true
    ic1=grad_slopes(start_condition,keep_grads);
    ic2=grad_slope_deltas(start_condition,keep_grads);
else % excite specific eigenmode
    mode_num=1;
    mode_column=(mode_num*2)-1;
    ic1=real(v(1:2:end,mode_column))';
    ic2=real(v(2:2:end,mode_column))';
end
ics=[];
for i=1:n_grads
    ics=[ics;ic1(i)];
    ics=[ics;ic2(i)];
end

% simulate gradient timeseries using ode45
% previous method i used, for comparison
if false
    [yp yv]=gradient_ode(ic1,ic2,coupling_params,pred_interval);
end

% position and velocity equations for n gradients
% based on real + imaginary parts of eigenvalues/eigenvectors
% rows are: gradient 1 position, gradient 1 velocity, ...
% columns are: eigenvalue 1 real, eigenvalue 1 imaginary, ...
% gradient position/velocity are a summation of the eigenvectors at each
% timepoint (real + imaginary parts), scaled by their eigenvalues
% eigenvalue real part = damping (exponential growth or decay of signal amplitude)
% eigenvalue imaginary part = frequency
% eigenvector weights = weight of each gradient position/velocity on each eigenvector
t=0:step_size:pred_interval;
no_damping=false; % zero out the damping effect (changes in amplitude over time)
pos_vel_mat=zeros(n_grads*2,n_grads*4,length(t));
eigvec_mat=zeros(n_grads*2,n_grads*4);
formula_strs={};
ks_keep=[];
for i=1:(n_grads*2) % one row at a time: each eigvector component
    vec_num=i;
    for j=1:2:(n_grads*4) % two columns at a time: real and imaginary for each eigenvalue
        val_num=((j+1)/2);
        cur_eigval=d(val_num,val_num);
        cur_eigvalreal=real(cur_eigval);
        if no_damping
            cur_eigvalreal=0;
        end
        cur_eigvalimag=imag(cur_eigval);
        cur_exp=exp(cur_eigvalreal.*t);
        cur_eigvec=v(vec_num,val_num);
        cur_eigvec_real=real(cur_eigvec);
        cur_eigvec_imag=imag(cur_eigvec);
        real_comp1=cur_eigvec_real*cos(cur_eigvalimag.*t);
        real_comp2=-cur_eigvec_imag*sin(cur_eigvalimag.*t);
        imag_comp1=cur_eigvec_imag*cos(cur_eigvalimag.*t);
        imag_comp2=cur_eigvec_real*sin(cur_eigvalimag.*t);
        pos_vel_mat(i,j,:)=cur_exp.*(real_comp1+real_comp2);
        pos_vel_mat(i,j+1,:)=cur_exp.*(imag_comp1+imag_comp2);
        if mod((j+1)/2,2) % only do for odd eigenvalues, first half of conjugate pair
            cur_formula_str=sprintf('e^%1.4ft * \n(k%d[%2.4fcos(%2.4ft) - %2.4fsin(%2.4ft)] + \n k%d[%2.4fcos(%2.4ft) + %2.4fsin(%2.4ft)])',...
                cur_eigvalreal,j,cur_eigvec_real,cur_eigvalimag,cur_eigvec_imag,cur_eigvalimag,...
                j+1,cur_eigvec_imag,cur_eigvalimag,cur_eigvec_real,cur_eigvalimag);
            if j>=3
                cur_formula_str=sprintf('%s%s%s',formula_strs{i},' +\n ',cur_formula_str);
            end
            formula_strs{i}=cur_formula_str;
        end
    end
end

% solve for k parameters given initial condtions at time t=0
% keep first half of each conjugate pair
ks_keep=[];
for i=2:2:(n_grads*4)
    if mod(i/2,2)
        ks_keep=[ks_keep;i-1;i];
    end
end
ks=linsolve(pos_vel_mat(:,ks_keep,1),ics);

% derive position/velocity timeseries by scaling solutions by k constants given initial conditions
clear pos_vel
for i=1:length(t)
    pos_vel(:,i)=pos_vel_mat(:,ks_keep,i)*ks;
end

% plot position/velocity timeseries
if false
    %plotf(pos_vel(1:2:end,:)')
    figure;hold on;
    for i=1:length(keep_grads)
        cur_ind=(i*2)-1;
        h1=plot(pos_vel(cur_ind,:)','Color',xkcd_colors(keep_grads(i),:),'Linewidth',4);
        hold on
    end
end

% measure frequency information
ang_freqs=diag(imag(d(1:2:end,1:2:end))); % number of 2pi radians covered in one data timestep
nat_periods=(2*pi)./ang_freqs; % number of timesteps to cover one full cycle
timestep_length=2; % length in seconds of one data timestep (TR)
freqs=1./(nat_periods*timestep_length); % number of cycles per second (Hz)
freqs_cycle=1./(nat_periods);
beat_period=1/((freqs(2)-freqs(1))/2);

% power spectrum
[p,f] = pspectrum(pos_vel(1:2:end,:)');
if false
    figure
    for i=1:length(keep_grads)
        plot(f/pi,abs(p(:,i)),'Linewidth',2)
        hold on
    end
end

% eigenvector magnitude represents amplitude of each gradient's fluctuation
% on a given eigenmode
% k parameters scale each eigenmode based on the initial conditions
% eigenvector magnitude * k = final amplitude for each gradient on each eigenmode
clear grad_amps
mags_real=real(v(:,1:2:end));
mags_imag=imag(v(:,1:2:end));
ks_ones=ones(12,1)*2;
for i=1:n_grads*2
    for j=1:n_grads
        cur_real=mags_real(i,j);
        cur_imag=mags_imag(i,j);
        cur_k_real=ks_ones((j*2)-1);
        cur_k_imag=ks_ones(j*2);
        grad_amps(i,j)=sqrt(((cur_k_real*(-cur_real+cur_imag))^2) + ((cur_k_imag*(cur_real+cur_imag))^2));
    end
end
%grad_pair_amps=grad_amps(1:2:end,:)*grad_amps(1:2:end,:)';

% derive fc matrix as sum of gradient covariance matrices within + across modes
% mode-wide components for each gradient
% indiv_mode_pos_vel dimensions:
% n_grads * 2 (position/velocity) X n_timepoints X n_modes
% (each mode in isolation with others zeroed out)
clear cur_mode_pos_vel indiv_mode_pos_vel
n_modes=n_grads;
for i=1:n_modes
    cur_mode_real_ind=(i*2)-1;
    cur_mode_imag_ind=cur_mode_real_ind+1;
    % keep ks for current mode; zero out the remaining ks
    cur_ks=zeros(length(ks),1);
    cur_ks(cur_mode_real_ind)=ks(cur_mode_real_ind);
    cur_ks(cur_mode_imag_ind)=ks(cur_mode_imag_ind);
    for j=1:length(t)
        cur_mode_pos_vel(:,j)=pos_vel_mat(:,ks_keep,j)*cur_ks;
    end
    indiv_mode_pos_vel(:,:,i)=cur_mode_pos_vel;
end

% pairwise covariance of gradients across all modes
% each gradient pair's covariance is sum of covariance of mode-specific
% components for each gradient
clear all_mode_grad_covs all_mode_grad_covs_mats
for i=1:n_grads
    % get all components for gradient i
    X1=squeeze(indiv_mode_pos_vel((i*2)-1,:,:));
    for j=1:n_grads
        % get all components for gradient j
        X2=squeeze(indiv_mode_pos_vel((j*2)-1,:,:));
        cur_cov=(X1'*X2)./(size(X1,1)-1);
        cur_cov_sum=sum(cur_cov(:));
        all_mode_grad_covs(i,j)=cur_cov_sum;
        all_mode_grad_covs_mats(:,:,i,j)=cur_cov;
    end
end

% verify that the following are identical:
% 1) fc matrix from simulated gradient covariance
% 2) fc matrix from sum of eigenmode gradient angles
grad_cov_1=all_mode_grad_covs;
grad_cov_2=cov(pos_vel(1:2:end,:)'); % gradient covariance from data summed across modes (eg overall gradient timeseries)
fcmat_1=zeros(246);
fcmat_2=zeros(246);
for i=1:n_grads
    for j=1:n_grads
        fcmat_1=fcmat_1+((grad_weights(1:246,i).*grad_weights(1:246,j)')*grad_cov_1(i,j));
        fcmat_2=fcmat_2+((grad_weights(1:246,i).*grad_weights(1:246,j)')*grad_cov_2(i,j));
    end
end
if false
    fcmat=fcmat_1;
    %plot_matrix;
    % confirmation that pairwise approach matches gradient covariance
    fcmat=fcmat_2;
    plot_matrix;
    k=find(triu(ones(246),1));
    plotcorrf(fcmat_1(k),fcmat_2(k))
end
fcmat_sim=fcmat_1;


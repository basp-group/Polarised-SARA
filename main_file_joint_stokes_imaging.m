% *************************************************************************
% Matlab code for sparse Stokes imaging under polarization constraint
% *************************************************************************
% 
% Related paper:
% Sparse interferometric Stokes imaging under polarization constraint 
% (Polarized SARA), MNRAS, 2018.
% Jasleen Birdi, Audrey Repetti, Yves Wiaux
%
% Contact: jb36@hw.ac.uk
% *************************************************************************

close all;
clear all;

%%%% Add folders to path
addpath('data');
addpath('lib');
addpath('Tools');
addpath('Algorithms');
addpath('TV');
addpath(genpath('irt'));


% -------------------------------------------------------------------------
% Parameters for the input image
% -------------------------------------------------------------------------

param_sim_data.sig_noise = 5e-3;
param_sim_data.Nx = 100;
param_sim_data.Ny = 100;
param_sim_data.N = param_sim_data.Nx*param_sim_data.Ny;
param_sim_data.im_choice = 'M87'; % file or M87
param_sim_data.im_type = 'counter_jet'; % forward_jet or counter_jet
param_sim_data.pixel_size = 14;
param_sim_data.im_fits = 0; % 1: image stored as fits file

% -------------------------------------------------------------------------
% sparsity regularisation
% -------------------------------------------------------------------------

param.dict = 'SARA'; %SARA, TV, Dirac+TV

% -------------------------------------------------------------------------
% sampling pattern parameters
% -------------------------------------------------------------------------
% options 'gaussian', 'file', 'gaussian+large-holes', 'file+undersample'
   
sampling_pattern = 'file';
param_sampling.file_name = 'data/avery_final.mat';
num_tests = 1; % number of tests to perform on the same coverage with different noise
% realisations


% Generate data
script_get_input_data_polar_RIME;


%%
% -------------------------------------------------------------------------
% Parameters for the Primal-dual algorithm
% -------------------------------------------------------------------------

% method
% 1: data fidelity + positivity
% 2: data fidelity + positivity + sparsity
% 3: data fidelity + positivity + sparsity + polarization constraint


if strcmp(param.dict,'SARA')
    proj_l2 = 1; % 1: Constrained formulation; 0: Unconstrained formulation
else
    proj_l2 = 0;
end

param_algo.tmax = 1e+5; % Maximum number of iterations
param_algo.Eta1 = 500;
param_algo.Eta2 = 2260;
param_algo.Eta3 = 1;
param_algo.tol_nnls = 1e-3; % Tolerance for l2 ball
param_algo.tol_rw = 5e-3;
param_algo.num_rw = 10; % Total number of reweights to be performed
param_algo.beta = 0.1; 
param_algo.sep = 1; % 0: projection onto single l2 ball; 1: projection onto separate balls
param_algo.bright = 1; % 0: Stokes matrix; 1:Brightness matrix
if proj_l2 == 1
method = 1; % For initialization to constrained formulation
param_algo.method_change = 3; % 2: without pol. const or  3: with pol. const.
param_algo.rel_stop_crit = 5e-5; % Relative variation between the solutions used for stopping criterion
param_algo.rel_stop_crit_change = 1e-5;
else 
    method = 3;
     param_algo.rel_stop_crit = 1e-5;
     if strcmp(param_sim_data.im_type,'counter_jet')
             param_algo.rel_stop_crit = 7e-6;
     end
end
rw = 0;
param_algo.tmax_rw = 10000;

param.pos = 0;
param.real = 0;
if strcmp(param_sim_data.im_type,'counter_jet')
    eta_o= [1e-5,1e-5,1e-5];  % Regularization parameter with l1 term
else
    eta_o = [1e-4, 1e-4, 1e-4];
end

parameters_algo;
param_algo.eta = eta_o;
param.gamma = (param_algo.eta(1)./sigma1);
param.verbose = 2; % 0: no log/images display;  1: to display log every 1000 iteration; 2: to display SNR plots, images every 1000 iterations; 
disp_cons = 0;


 
%%
% *****************************************************************************************
% Compute the solution
% Run primal-dual algorithm
% *****************************************************************************************

fprintf('Sparsity dictionary: %s\n\n', param.dict);
fprintf('Method: %d\n\n', method);

for n_test = 1:num_tests
def_var; % Define and initialize variables

fprintf('Starting algorithm:\n\n');

tstart = tic;
if bright
    pdfb_algo_stokes_imaging_bright; 
else
    pdfb_algo_stokes_imaging; 
end
tend = toc(tstart);

fprintf('Algorithm runtime: %ds\n\n', ceil(tend));
end





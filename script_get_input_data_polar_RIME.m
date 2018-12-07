% script that generates the input data for the test
% set random nr gen
rng shuffle;

    fprintf('Generating new data ... \n\n');
    %% image loading
    param_sim_data.P = 3; % Number of Stokes images considered   

    [im, N, Ny, Nx] = gen_image(param_sim_data);
    S = [im{1}(:), im{2}(:), im{3}(:)]; % Stokes matrix
    L = [1,0,0,1; 1,0,0,-1; 0,1,1,0]; % Conversion matrix
    Lt = 0.5*(conj(L))'; % Adjoint conversion matrix
    
    B = cell(4,1);
    B1 = S*L; % Brightness matrix
    
    
    for i = 1:4
        B{i} = reshape(B1(:,i), Ny, Nx);
    end
    
  
    %% generate the sampling pattern
    
    ox = 2; % oversampling factors for nufft
    oy = 2; % oversampling factors for nufft
    Kx = 8; % number of neighbours for nufft
    Ky = 8; % number of neighbours for nufft

    % for randomly generated coverages
    param_sampling.hole_number = 8000; % number of holes to introduce for 'gaussian+large-holes'
    param_sampling.hole_size = pi/60; % size of the missing frequency data
    param_sampling.sigma = pi/4; % variance of the gaussion over continous frequency
    param_sampling.sigma_holes = pi/3; % variance of the gaussion for the holes
    
    % data splitting config
    param_sampling.equal_partitioning = 1; % flag
    param_sampling.equal_partitioning_no = 1;
    
    param_sampling.N = param_sim_data.N; % number of pixels in the image
    param_sampling.Nox = ox*param_sim_data.Nx; % number of pixels in the oversampled image
    param_sampling.Noy = oy*param_sim_data.Ny; % number of pixels in the oversampled image
    param_sampling.pixel_size = param_sim_data.pixel_size;
    
    [uw, vw, ~] = util_gen_sampling_pattern(sampling_pattern, param_sampling);
    
    figure, plot(uw{1,1},vw{1,1},'.')
    
    u = uw;
    v = vw;
    
    nW = ones(length(uw), 1); % weights
     
    %% measurement operator initialization
    
    fprintf('Initializing the NUFFT operator\n\n');
    tstart = tic;
    [A, At, G, W, Gw] = op_p_nufft([v u], [Ny Nx], [Ky Kx], [oy*Ny ox*Nx], [Ny/2 Nx/2]);
    tend = toc(tstart);
    
    fprintf('Initialization runtime: %ds\n\n', ceil(tend));
    R = length(v);
    
    y0 = cell(num_tests, 1);
    yb = cell(num_tests, 1);
    sigma_noise = cell(num_tests, 1);
    
    %% generate noisy input data
  
    
    for seed = 1:num_tests
        
        n_corr = 4; % Number of considered correlations
 
        [y0{seed}, yb{seed}, sigma_noise{seed}] = util_gen_input_data_RIME(B, Gw, A, param_sim_data.sig_noise,3,n_corr,seed);
        
        y_conc = [yb{seed,1}{1:n_corr,1}];
        y_stokes = y_conc*Lt;
        ys{seed} = (mat2cell(y_stokes, size(y_stokes,1), ones(1,size(y_stokes,2))))'; % Unmixed measurement vector
        
    end
    
    y = ys;
    yT = y;

    %% sparsity operator definition
    
    nlevel = 4; % wavelet level
    
    if strcmp(param.dict,'TV')
        Psitw = @(x) ComputeTVlin(reshape(x,Nx,Ny));
        Psiw = @(x) ComputeTVlinAdj(x,Nx,Ny);
    elseif strcmp(param.dict,'Dirac+TV')
        Psit_tv = @(x) ComputeTVlin(reshape(x,Nx,Ny));
        Psi_tv = @(x) ComputeTVlinAdj(x,Nx,Ny);
        
        [Psiw, Psitw] = SARA_sparse_operator(im{1}, nlevel,'Dirac');
    else
        [Psiw, Psitw] = SARA_sparse_operator(im{1}, nlevel,param.dict);
    end
    
    param_sim_data.Psit = Psitw;
    param_sim_data.Psi = Psiw;
    
    P = param_sim_data.P;
    Psitw = cell(1,P);
    Psiw = cell(1,P);
    for i=1:P
        Psitw{i} = param_sim_data.Psit;
        Psiw{i} =  param_sim_data.Psi;
    end
    
    %% compute the operator norm
    fprintf('Computing operator norm ...\n');
    
    fprintf('Natural W ...\n');
    verbosity = 0;
    evl = op_norm(@(x) Gw * A(x), @(x) At(Gw' * x), [Ny, Nx], 1e-4, 200, verbosity);

    %% L2 ball sizes
    epsilon = cell(num_tests, 1);
    epsilons = cell(num_tests, 1);
    epsilonT = cell(num_tests, 1);
    epsilonTs = cell(num_tests, 1);
    
    
    %% generate bounds
    %% Projection onto l2 ball

% options: 
% l2_ball_definition -> 'sigma', 'chi-percentile', 'value'
% stopping_criterion -> 'sigma', 'chi-percentile', 'l2-ball-percentage', 'value'

l2_ball_definition = 'sigma';
stopping_criterion = 'l2-ball-percentage';

param_l2_ball.stop_eps_v = [sqrt(2*307780)];
param_l2_ball.val_eps_v = 1.04 * param_l2_ball.stop_eps_v;

param_l2_ball.sigma_ball = 2;
param_l2_ball.sigma_stop = 2;

param_l2_ball.chi_percentile_ball = 0.99;
param_l2_ball.chi_percentile_stop = 0.999;

param_l2_ball.l2_ball_percentage_stop = 1.0001;

    for k = 1:num_tests
%         
            [eps{k}] = util_gen_L2_bounds_polar(yb{k}, param_sim_data.P, ...
                param_sim_data.sig_noise);
           
            epsilon{k} = eps{1,1}(1); % Noise bounds for each correlator
    end
    
    

% Define variables for the algorithm
P = param_sim_data.P;

St_im = cell(P+2,1); % Stokes images: I,Q,U and constraint variables zeta

% dummy dual variables
dual_tilde1 = cell(P,1); 
dual_tilde2 = cell(P+2,1); 
dual_tilde3 = cell(P,1);
dual_tilde4 = cell(P,1);

% dual variables
dual_var1 = cell(P,1); 
dual_var2 = cell(P+2,1); 
dual_var3 = cell(P,1);
dual_var4 = cell(P,1);

dummy1 = cell(P+1,1);
dummy2 = cell(P+1,1);
dummy3 = cell(P+1,1);
dummy4 = cell(P+1,1);

St_im_old = cell(P+2,1);
epi_h = cell(2,1);

snr_v = cell(P,1);
mse = cell(P,1);
nrmse = cell(P,1);
rel_sol = cell(P,1);
l1_norm = zeros(param_algo.tmax,1);
residual = zeros(param_algo.tmax,1);
St_iter = cell(P,param_algo.tmax);

N = Nx*Ny;

sep = param_algo.sep;
bright = param_algo.bright;
dict = param.dict;
if proj_l2 == 1
    tol = param_algo.tol_nnls;
else
    tol = param_algo.tol_rw;
end

if bright
    dual_var3 = cell(4,1);
    if sep ~= 1
        epsilon{n_test} = epsilon{n_test}*sqrt(4);
    end
else
    epsilon{n_test} = epsilon{n_test}/sqrt(2); % Noise bound for Stokes visibility
end

init_var; % To initialize the variables

t_prev = 0;


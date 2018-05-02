% Initialize variables

for i = 1:P+2

    
  St_im{i} = zeros(Nx,Ny);

 dual_tilde1{i} = zeros(size(param_sim_data.Psit(St_im{i})));
 dual_var1{i} = zeros(size(param_sim_data.Psit(St_im{i})));

 dual_tilde2{i} = zeros(size(St_im{i}));
 dual_var2{i} = zeros(size(St_im{i}));
 
 dual_tilde3{i} = zeros(size(Gw*A(St_im{i})));
 dual_var3{i} = zeros(size(Gw*A(St_im{i})));
 
 if strcmp(dict,'Dirac+TV')
     dual_tilde4{i} = zeros(size(St_im{i}));
     dual_var4{i} = zeros(size(Psit_tv(St_im{i})));
 end

%  snr_v{i} = zeros(tmax, 1);
%  mse{i} = zeros(tmax, 1);
%  nrmse{i} = zeros(tmax, 1);
 rel_sol{i} = zeros(param_algo.tmax,1);
 
 
end
% 

for i = 1:P
    g1{i} = zeros(size(St_im{i}));
    g2{i} = zeros(size(St_im{i}));
    g3{i} = zeros(size(St_im{i}));
end

if method >1
epi_h{1} = 0*[dual_var2{1}(:),dual_var2{4}(:)];
epi_h{2} = 0*[dual_var2{2}(:),dual_var2{3}(:), dual_var2{5}(:)];
end

if bright && sep~=1
dual_tilde3 = zeros(length(y{1,n_test}{i,1}),4);
dual_var3 = zeros(length(y{1,n_test}{i,1}),4);
end


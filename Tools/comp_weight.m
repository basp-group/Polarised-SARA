% Script to compute weights

temp = param_sim_data.Psit(reshape(St_im{i},Nx,Ny));

sigma_n = param_sim_data.sig_noise*sqrt(numel(y{1,n_test}{i,1})/(N*9));

% Weights
if strcmp(dict,'TV')
    u = temp(1:Nx*Ny);
    v = temp(Nx*Ny+1:end);
   temp = sqrt(u.^2 + v.^2);
end
    weights{i}=abs(temp);

sigma_s=std(temp(:));

if rw == 1
  delta{i}(rw) = max(temp(:));
end

if rw >1
    delta{i}(rw) = delta{i}(rw-1)*param_algo.beta;
end


weights{i}=delta{i}(rw)./(delta{i}(rw)+weights{i});

    
    
    

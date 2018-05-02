function [y0, y, sigma_noise] = util_gen_input_data_RIME(B,Gw, A, sigma_noise,P,pp,seed)
% generates the input data


y0 = cell(pp,1);
y = cell(pp,1);
% sigma_noise = cell(pp,1);
% normy0 = cell(pp,1);

rng(seed)

for i = 1:pp
    
y0{i} = Gw * A(B{i}); 
Nm = numel(y0{i});
% normy0{i} = norm(y0{i});


% add Gaussian noise
% sigma_noise{i} = 10^(-input_snr(i)/20) * normy0{i}/sqrt(Nm);

noise = (randn(Nm, 1) + 1i*randn(Nm, 1)) * sigma_noise/sqrt(2);
y{i} = y0{i} + noise;
end


end

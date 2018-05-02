function [y0, y, sigma_noise] = util_gen_input_data(im, Gw, A, input_snr,P)
% generates the input data


y0 = cell(P,1);
y = cell(P,1);
sigma_noise = cell(P,1);
normy0 = cell(P,1);


for i = 1:P
    
y0f = A(im{i});
y0{i} = Gw * y0f; 
Nm = numel(y0{i});
normy0{i} = norm(y0{i});

% add Gaussian i.i.d. noise
sigma_noise{i} = 10^(-input_snr(i)/20) * normy0{i}/sqrt(Nm);

noise = (randn(Nm, 1) + 1i*randn(Nm, 1)) * sigma_noise{i}/sqrt(2);
y{i} = y0{i};
y{i} = y{i} + noise;
end


end


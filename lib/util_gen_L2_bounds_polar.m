function [epsilon] = util_gen_L2_bounds_polar(y,P,sigma_noise)
% generates the l2 noise bounds

Nm = numel(y{1,1});

for i = 1:P
 
y0 = y{i,1};  

R = length(y{i,1});
Nm = numel(y{i,1});
% normy0 = norm(y{i,1});
epsilon(i) = sqrt(2*Nm + 2*sqrt(2*2*Nm)) * sigma_noise/sqrt(2);

end

end


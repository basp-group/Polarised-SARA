function [V] = Proj_epih1(sigma,v)
% Epigraphical projection on h1

v = v./sigma;

v1 = v(:,1);
v2 = v(:,2);

h = ones(length(v1),1);
h(v1+v2<0) = 0;

v1(~h) = 0.5*(v1(~h) - v2(~h));
% v2(~h) = max(0.5*(v2(~h)-v1(~h)), v2(~h));

v2(~h) = max(-v1(~h), v2(~h));

V = [v1,v2];

% To check the constraint
cc = size(find(-v1>v2));

if cc(1) > 0
    fprintf('Dual_Epi_h1 not satisfied:%d',cc(1));
end


end
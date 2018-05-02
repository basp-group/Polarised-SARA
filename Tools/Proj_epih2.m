function [V] = Proj_epih2(sigma,v)
% Epigraphical projection on h2

v = v./sigma;

v1 = v(:,1);
v2 = v(:,2);
v3 = v(:,3);

V = [v1,v2,v3];

v12_norm = sqrt(abs(v1).^2 + abs(v2).^2);

h = ones(length(v12_norm),1);
h(v12_norm<-v3) = 0;

V(~h,:) = 0;

h(v12_norm < v3) = 0;

v_frac = v3./v12_norm;
v_frac(isnan(v_frac)) = 0;
alpha = 0.5*(ones(length(v1),1) + v_frac);
% V(h,:) = alpha.*[v1,v2,v12_norm];
v_dummy = bsxfun(@times,[v1,v2,v12_norm],alpha);

h(h==0) = 2;
h(h==1) = 0;
% V(~h2,:) = v_dummy(~h2,:);
V(~h,:) = v_dummy(~h,:);

% To check the constraint

tol_r = 1e-8;

v12 = sqrt(abs(V(:,1)).^2 + abs(V(:,2)).^2);
r = V(:,3)./v12;

cc = size(find(r < 1-tol_r));

if cc(1) > 0
    fprintf('Dual_Epi_h2 not satisfied:%d \n',cc(1));
end
end

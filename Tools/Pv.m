
function [Z] = Pv(z)

% Projection on the constraint set
z1 = z(:,1);
z2 = z(:,2);

Z = z;

z_idx = ones(length(z1),1); % indices for constraint set 
z_idx(z1+z2>0) = 0;
Z(~z_idx,:) = bsxfun(@minus,z(~z_idx,:),(z1(~z_idx)+z2(~z_idx)));

end

    
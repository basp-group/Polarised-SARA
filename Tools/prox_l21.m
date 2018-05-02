function [B_res] = prox_l21(B,lambda)

% B is the matrix whose l2,1 norm needs to be computed
% lambda is the associated thresholding parameter

N = size(B,1);
% B = [real(B), imag(B)];
% a1 = B(1:end,1:n);
% a2 = B(1:end,n+1:end);
% B = [B(:,1), B(:,2)];
B_res = zeros(size(B)); 


for n = 1:N
    if norm(B(n,:)) >= lambda 
  B_res(n,:) =  B(n,:).*((norm(B(n,:)) - lambda)./norm(B(n,:)));
    end
end

% B_res = B_res(:,1) + 1i*B_res(:,2);
% B_res = [reshape(B_res(:,1),n,n), reshape(B_res(:,2),n,n)];
B_res = [B_res(:,1), B_res(:,2)];


function [res_norm, DR, DR_peak_rms] = comp_dynamic_range(N, A, Gw, At, St_im, evl, y, p)
% Compute the dynamic range

x = St_im{p};
peak(p) = max(x(:));
rms_n(p) = rms(St_im{p}(:));
num(p) = sqrt(N)* evl * max(x(:));

den_term = real(At(Gw'*((Gw*A(x)-y{p,1}))));
res{p} = den_term;
den(p) = norm(den_term,2);
res_norm(p) = den(p)./(sqrt(N)* evl); % Estimated noise level in image x
DR(p) = num(p)./den(p);

DR_peak_rms(p) = peak(p)./rms_n(p);

end




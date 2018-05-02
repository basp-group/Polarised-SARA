
% Compute the pixels not satisfying polarization constraint, lying over
% noise threshold

pol_err = P_abs - St_im{1};
mask = zeros(size(St_im{1}));
mask(pol_err>0) = pol_err(pol_err>0); % Map of non-physical valued pixels before thresholding
p = 1; % noise to be estimated from total intensity image
[res_norm, ~, ~] = comp_dynamic_range(param_sim_data.N, A, Gw, At, St_im, evl, y{1,n_test}, p);

% dynamic_range_calc; % Perform thresholding

mask2 = zeros(size(St_im{1}));
mask2(pol_err>(3*res_norm(1))) = pol_err(pol_err>(3*res_norm(1))); %Map of non-physical valued pixels after thresholding

pol_thresh =  size(find(mask2(:)),1); % Number of non-physical valued pixels

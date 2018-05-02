function [mse,nrmse] = calculate_mse(im,St_im)
% Computes MSE and NRMSE
% im: True image
% St_im: Reconstructed image

err = abs(St_im-im).^2;
err = sum(err(:));
n = abs(im).^2;
n = sum(n(:));
mse = err;
nrmse = sqrt(mse/n);

end
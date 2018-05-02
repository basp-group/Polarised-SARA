% Generate image

im_ex = xlsread('avery_sgra_excel_file');

im = im_ex(9:end,3); % Stokes I
% imagesc(im)
s = size(im,1);
im = reshape(im,sqrt(s),sqrt(s));
% imagesc(im)

if resize
im = padarray(im,[14,14]);
end
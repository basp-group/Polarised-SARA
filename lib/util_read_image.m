function [im, N, Ny, Nx] = util_read_image(name, normalize)
% Reads and takes care of retriving all the info from an fits image
%
% in:
% name   - image name
%
% out:
% im     - matrix representing the image
% N      - number of pixels Ny*Nx
% Ny, Ny - image imensions

fits_file = 1;


    
if ~exist('normalize', 'var')
    normalize = 1;
end


if fits_file == 1
data_raw = fitsread(name);
data = flipud(data_raw);

else
    data = name;
end


im = im2double(data);

[Ny, Nx] = size(im); 
N = Ny*Nx;

if normalize == 1
    im(im<0) = 0;
    % scale image to avoid numerical errors
    im = im/max(max(im));
end

end


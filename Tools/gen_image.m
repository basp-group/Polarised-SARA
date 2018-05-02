function [im, N, Ny, Nx] = gen_image(param_sim_data)
% Reads image
%
% in:
% param_sim_data  - structure containing desired image specifications
%
% out:
% im     - cell representing the Stokes images
% N      - number of pixels Ny*Nx
% Ny, Nx - image dimensions

if strcmp(param_sim_data.im_choice,'file')
    
    im_ex = xlsread('avery_sgra_excel_file');
    im = cell(3,1);
    
    for i = 1:param_sim_data.P
        im{i} = im_ex(9:end,i+2);
        s = size(im{i},1);
        im{i} = reshape(im{i},sqrt(s),sqrt(s));
    end
    
    im{1}(im{1}<0) = 0;
    [Ny, Nx] = size(im{1});
    N = Ny*Nx;
    
elseif ~param_sim_data.im_fits
    
    Nx = param_sim_data.Nx;
    Ny = param_sim_data.Ny;
    N = Ny*Nx;
    
    if strcmp(param_sim_data.im_type,'forward_jet')
        ii = fitsread('avery_hm2.I.model.fits');
        qq = fitsread('avery_hm2.Q.model.fits');
        uu = fitsread('avery_hm2.U.model.fits');
        
        im{1} = ii(230:295,220:305);
        im{2} = qq(230:295,220:305);
        im{3} = uu(230:295,220:305);
        
    else
        ii = fitsread('jason-j2.I.model.fits');
        qq = fitsread('jason-j2.Q.model.fits');
        uu = fitsread('jason-j2.U.model.fits');
        
        im{1} = ii(60:240,60:240);
        im{2} = qq(60:240,60:240);
        im{3} = uu(60:240,60:240);
        
    end
    
    im{1} = imresize(im{1},[Nx,Ny]);
    im{2} = imresize(im{2},[Nx,Ny]);
    im{3} = imresize(im{3},[Nx,Ny]);
    
    im{1}(im{1}<0) = -im{1}(im{1}<0);
    
    
else
    for i = 1:param_sim_data.P
        im{i} = flipud(fitsread(param_sim_data.im_choice));
        im{i} = im2double(im{i});
    end
    
    im{1}(im{1}<0) = 0;
    [Ny, Nx] = size(im{1});
    N = Ny*Nx;
end

end





% Plot images showing pixels not satisfying flux constraint

method = 1;
seed = 2;

% Choice of dictionary
% 1: 'TV'
% 2: 'Dirac+TV'
% 3: 'SARA'

flux_err1_avery = zeros(3,1);
flux_err1_thresh_avery = zeros(3,1);
flux_err3_avery = zeros(3,1);
flux_err3_thresh_avery = zeros(3,1);
flux_err1_jason = zeros(3,1);
flux_err1_thresh_jason = zeros(3,1);
flux_err3_jason = zeros(3,1);
flux_err3_thresh_jason = zeros(3,1);


for k = 1:3
    switch k
        case 1
            dict = 'TV';
        case 2
            dict = 'Dirac+TV';
        case 3
            dict = 'SARA';
    end
    
    for method = 1:2:3
    name_file = sprintf('Save_rec_dict=%s,type=jason,method=%d,seed=%d.mat',dict,method,seed);
    
    load(name_file)




% Define variables
% I_true = im{1};
% Q_true = im{2};
% U_true = im{3};
% P_true = Q_true + 1i* U_true;
% I_rec = St_im{1};
% Q_rec = St_im{2};
% U_rec = St_im{3};
% P_rec = Q_rec + 1i*U_rec;

I_true = im{1};
Q_true = im{2};
U_true = im{3};
P_true = Q_true + 1i* U_true;
% I_rec = St_im_b{1,2};
% Q_rec = St_im_b{2,2};
% U_rec = St_im_b{3,2};
I_rec = St_im_f{1};
Q_rec = St_im_f{2};
U_rec = St_im_f{3};
P_rec = Q_rec + 1i*U_rec;

% 
% 
%%

% Compute flux constraint error images
P_abs = sqrt(abs(Q_rec).^2+abs(U_rec).^2);
flux_err = P_abs - I_rec;

mask = zeros(size(I_true));
% mask(flux_err>0) = flux_err(flux_err>0);
% figure, imagesc(mask)
mask(flux_err>0) = flux_err(flux_err>0);
figure,imagesc(mask)
colorbar('FontSize',16,'FontWeight','bold')
axis image, axis off

if method == 1
flux_err1_jason(k) = size(find(mask(:)),1);
fprintf('Number of pixel with unsat. const.= %d\n', flux_err1_jason(k));
else
   flux_err3_jason(k) = size(find(mask(:)),1);
fprintf('Number of pixel with unsat. const.= %d\n', flux_err3_jason(k));
end

% 
% figure
% subplot 221, imagesc(St_im{1}), axis image, axis off
% subplot 222, imagesc(P_abs), axis image, axis off
% subplot 223, imagesc(mask), axis image, axis off
% subplot 224, imagesc(log10(P_abs)), axis image, axis off


%% Perform thresholding

dynamic_range_calc;

mask2 = zeros(size(I_true));
mask2(flux_err>(3*res_norm(1))) = flux_err(flux_err>(3*res_norm(1)));
figure, imagesc(mask2)
colorbar('FontSize',16,'FontWeight','bold')
axis image, axis off

if method == 1
flux_err1_thresh_jason(k) = size(find(mask2(:)),1);
fprintf('Number of pixel with unsat. const. after threshold= %d\n', flux_err1_thresh_jason(k));
else
   flux_err3_thresh_jason(k) = size(find(mask2(:)),1);
fprintf('Number of pixel with unsat. const. after threshold= %d\n', flux_err3_thresh_jason(k));
end 

%% RMS values 

% ri = rms(I_rec(:));
% rp = rms(P_rec(:));
% 
% % Polarization fraction images
% 
% % True images
% m = P_abs./I_rec;
% m_t = zeros(size(I_true));
% m_t(I_rec>2*ri) = m(I_rec>2*ri);
% figure
% subplot 121, imagesc(m), axis image
% subplot 122, imagesc(m_t), axis image
    end
end

save('flux_error_counts.mat','flux_err1_jason','flux_err3_jason',...
    'flux_err1_thresh_jason','flux_err3_thresh_jason')
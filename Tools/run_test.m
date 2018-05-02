% 

% Define variables
P = 3;

St_im = cell(P+2,1); % Stokes images: I,Q,U and constraint variables zeta
dual_tilde1 = cell(P,1); % dummy dual variables
dual_tilde2 = cell(P+2,1); % dummy dual variables

dual_var1 = cell(P,1); % dual variables
dual_var2 = cell(P+2,1); % dual variables

dummy1 = cell(P+1,1);
dummy2 = cell(P+1,1);
St_im_old = cell(P+2,1);
epi_h = cell(2,1);

snr_v = cell(P,1);
rel_sol = cell(P,1);
l1_norm = zeros(tmax,1);
residual = zeros(tmax,1);

N = Nx*Ny;

% Initialize variables


for i = 1:P+2

 St_im{i} = zeros(Nx,Ny);
%  St_im{i} = im{i};
 
 dual_tilde1{i} = zeros(size(Psitw(St_im{i})));
 dual_var1{i} = zeros(size(Psitw(St_im{i})));

 dual_tilde2{i} = zeros(size(St_im{i}));
 dual_var2{i} = zeros(size(St_im{i}));
 
 snr_v{i} = zeros(tmax, 1);
 rel_sol{i} = zeros(tmax,1);
 
end

epi_h{1} = 0*[dual_var2{1}(:),dual_var2{4}(:)];
epi_h{2} = 0*[dual_var2{2}(:),dual_var2{3}(:), dual_var2{5}(:)];

%  
%  if i == 1
%  zeta1{i} = 0;
%  zeta2{i} = 0;
% 
%  dual_zeta1_tilde{i} = 0;
%  dual_zeta2_tilde{i} = 0;
% 
%  dual_zeta1{i} = 0;
%  dual_zeta2{i} = 0;
% 
%  end



%% Functions for projection

% thresholding negative values
negt = @(z) max(real(z), 0);

%soft thresholding operator
soft = @(z, T) sign(z) .* max(abs(z)-T, 0); 

%% Algorithm

for t = 1:tmax
    % In parallel
    
    for i = 1:P+1
        St_im_old{i} = St_im{i};
        if i <= P
      dummy1{i} = St_im{i} - tau * real(At(Gw'*(Gw*A(St_im{i})-y{1,1}{i,1})));  
      
      
      dummy2{i} = tau*(Psiw(dual_var1{i}) + dual_var2{i});
      

      St_im{i} = dummy1{i} - dummy2{i};
      if i == 1
          St_im{i} = negt(St_im{i});
      end
      
        dual_tilde1{i} = dual_var1{i} + sigma1*Psitw(2*St_im{i}-St_im_old{i});
        dual_var1{i} = dual_tilde1{i} - sigma1*soft(dual_tilde1{i}./sigma1, eta(i)/sigma1);
       
        snr_v{i}(t) = 20*log10(norm(im{i}(:))/norm(im{i}(:) - St_im{i}(:)));
        else
            St_im_old{i+1} = St_im{i+1};
            dummy1{i} = [St_im{i}(:),St_im{i+1}(:)];
            dummy2{i} = tau.*[dual_var2{i}(:), dual_var2{i+1}(:)];
            St{i} = Pv(dummy1{i}-dummy2{i});
            St_im{i} = reshape(St{i}(:,1),Nx,Ny);
            St_im{i+1} = reshape(St{i}(:,2),Nx,Ny);

        end
        
        dual_tilde2{i} = dual_var2{i} + sigma2 * (2*St_im{i}-St_im_old{i});
        
    end
    
    dual_tilde2{P+2} = dual_var2{P+2} + sigma2 * (2*St_im{P+2}-St_im_old{P+2});

    
    for i = 1:2
       if i == 1
           vec = [dual_tilde2{1}(:),dual_tilde2{4}(:)];
           epi_h{i} = vec - sigma2*Proj_epih1(sigma2,vec);
       
       else
        vec = [dual_tilde2{2}(:),dual_tilde2{3}(:),dual_tilde2{5}(:)];
        epi_h{i} = vec - sigma2*Proj_epih2(sigma2,vec);
       end
    end
    
    % Update variables
    dual_var2{1} = reshape(epi_h{1}(:,1),Nx,Ny);
    dual_var2{4} = reshape(epi_h{1}(:,2),Nx,Ny);
    
    dual_var2{2} = reshape(epi_h{2}(:,1),Nx,Ny);
    dual_var2{3} = reshape(epi_h{2}(:,2),Nx,Ny);
    dual_var2{5} = reshape(epi_h{2}(:,3),Nx,Ny);

    rm = 0;
    for i = 1:P
        
    l1_norm(t) = l1_norm(t) + eta(i)*sum(abs(Psitw(St_im{i})));
    residual(t) = residual(t) + 0.5*sum((abs(Gw*A(St_im{i})-y{1,1}{i,1})).^2);
    rel_sol{i}(t) = norm(St_im{i}(:) - St_im_old{i}(:))/norm(St_im_old{i}(:));
    
    rm = max(rel_sol{i}(t),rm);
    end
    
    obj(t) = l1_norm(t) + residual(t) ;
    
    %% Stopping criterion and logs
    
    % log
    if (param.verbose >=4)
        fprintf('Iter %i\n',t);
        fprintf('L1 norm              = %e\n',l1_norm(t))
        fprintf('Residual             = %e\n',residual(t))
        fprintf('Obj function value   = %e\n',obj(t))
%         fprintf('Rel sol norm change  = %e\n',rel_sol)
    end

    if (param.verbose == 2)
        figure(100)
        subplot 331, semilogy(obj), title('Obj function')
        subplot 334, plot(snr_v{1}), title('SNR- I')
        subplot 335, plot(snr_v{2}), title('SNR- Q')
        subplot 336, plot(snr_v{3}), title('SNR- U')
        subplot 337, semilogy(rel_sol{1}), title('rel norm- I')
        subplot 338, semilogy(rel_sol{2}), title('rel norm- Q')
        subplot 339, semilogy(rel_sol{1}), title('rel norm- U')
        
        figure(101)
        subplot 231, imagesc(im{1}), title('True I'), axis image, axis off
        subplot 232, imagesc(im{2}), title('True Q'), axis image, axis off
        subplot 233, imagesc(im{3}), title('True U'), axis image, axis off
        subplot 234, imagesc(St_im{1}), title('Reconstructed I'), axis image, axis off
        subplot 235, imagesc(St_im{2}), title('Reconstructed Q'), axis image, axis off
        subplot 236, imagesc(St_im{3}), title('Reconstructed U'), axis image, axis off
        
    pause(0.1)
    end
    
    
       if rm <= 1e-4
           break
       end
end
    

% if reweighting
%     reweight_polar;
% end
% 
        SNR = (snr_v{1}(end)+snr_v{2}(end)+snr_v{3}(end))/3;


%% Functions used in proximity operators

% scaling, projection on L2 norm
sc = @(z, radius) z * min(radius/norm(z(:)), 1);

% thresholding negative values
negt = @(z) max(real(z), 0);

%soft thresholding operator
soft = @(z, T) sign(z) .* max(abs(z)-T, 0);

% Positive soft thresholding operator
pos_soft = @(z, T) max(real(z)-T, 0);

%%
%****************************************************************************************%
%****************************************************************************************%

% Algorithm

for t = 1:param_algo.tmax
    tm = tic;
    % par
    for i=1:P
        %% Primal update
        St_im_old{i} = St_im{i};
        
        if method == 1
            St_im{i} = real(St_im{i} - tau*(reshape(g3{i},Ny,Nx)));
        elseif method == 2
            St_im{i} = real(St_im{i} - tau*(reshape(g1{i},Ny,Nx) + reshape(g3{i},Ny,Nx)));
        else
            St_im{i} = real(St_im{i} - tau*(reshape(g1{i},Ny,Nx) + reshape(g2{i},Ny,Nx)+reshape(g3{i},Ny,Nx)));
        end
        if i == 1
            St_im{i} = negt(St_im{i});
        end
        
        s_b{i} = (2*St_im{i}-St_im_old{i});
        
        %% Dual update
        
        %% L1 prox update
        
        if method ~= 1
            
            dual_tilde1{i} = dual_var1{i} + sigma1*Psitw{i}(2*St_im{i}-St_im_old{i});
            
            if strcmp(dict,'TV')
                dual_var1{i} = dual_tilde1{i} - sigma1*prox_tv(1,dual_tilde1{i}./sigma1, eta(i)/sigma1,Nx,Ny);
                g1{i} = Psiw{i}(dual_var1{i});
                
            elseif strcmp(dict,'Dirac+TV')
                dual_tilde4{i} = dual_var4{i} + sigma4*Psit_tv(2*St_im{i}-St_im_old{i});
                dual_var1{i} = dual_tilde1{i} - sigma1*soft(dual_tilde1{i}./sigma1, eta(i)/sigma1);
                dual_var4{i} = dual_tilde4{i} - sigma4*prox_tv(1,dual_tilde4{i}./sigma4, nu(i)/sigma4,Nx,Ny);
                g1{i} = Psiw{i}(dual_var1{i}) + reshape(Psi_tv(dual_var4{i}),Nx,Ny);
            else
                dual_var1{i} = dual_tilde1{i} - sigma1*soft(dual_tilde1{i}./sigma1, eta{i}/sigma1);
                g1{i} = Psiw{i}(dual_var1{i});
            end
        end
        
         % Epi_proj
       dual_tilde2{i} = dual_var2{i} + sigma2 * (2*St_im{i}-St_im_old{i});
    end
    
    
    %% Data fidelity update
    
    s_m = [s_b{1}(:),s_b{2}(:), s_b{3}(:)];
    bb = s_m*L;
    b1 = cell(4,1);
    
    parfor j = 1:4
        Bb{j} = reshape(bb(:,j), Ny, Nx);
        
        phi_s{j} = Gw*A(Bb{j}); % To simulate the operation of phi(S)
        
        if sep == 1
            dual_tilde3{j} = dual_var3{j} + sigma3*phi_s{j};
            b1{j} = sc((dual_tilde3{j}./sigma3)-yb{n_test,1}{j,1},epsilon{n_test})+yb{n_test,1}{j,1};
            dual_var3{j} = dual_tilde3{j} - sigma3*b1{j};
            d3{j} = (At(Gw'*(dual_var3{j})));
            d3{j} = d3{j}(:);
        end
    end
    
    if sep~=1
        for j =1:4
            dual_tilde3(:,j) = dual_var3(:,j) + sigma3*phi_s{j};
            d3{j} = (At(Gw'*(dual_var3(:,j))));
            d3{j} = d3{j}(:);
        end
        
        dt3 = dual_tilde3(:);
        b1 = sc((dt3./sigma3) - y_conc(:), epsilon{n_test}) + y_conc(:);
        bb = vec2mat(b1',5631);
        bb = bb';
        dual_var3 = dual_tilde3 - sigma3*bb;
    end
    
    g3 = mat2cell([d3{1,1:4}]*Lt,[param_sim_data.N], [1 1 1])';
    %
    
    
    %% Epigraphical projection
    if method == 3 % Sparsity + polarization constraint
        i = 4;
        
        %Primal update
        St_im_old{i} = St_im{i};
        St_im_old{5} = St_im{5};
        g11{i} = [St_im{i}(:),St_im{i+1}(:)];
        g22{i} = tau*[dual_var2{i}(:), dual_var2{i+1}(:)];
        St{i} = Pv(g11{i}-g22{i});
        St_im{i} = reshape(St{i}(:,1),Ny,Nx);
        St_im{i+1} = reshape(St{i}(:,2),Ny,Nx);
        dual_tilde2{i} = dual_var2{i} + sigma2 * (2*St_im{i}-St_im_old{i});
        dual_tilde2{5} = dual_var2{5} + sigma2 * (2*St_im{5}-St_im_old{5});
        
        %Dual update
        vec = [dual_tilde2{1}(:),dual_tilde2{4}(:)];
        epi_h{1} = vec - sigma2*Proj_epih1(sigma2,vec);
        
        vec = [dual_tilde2{2}(:),dual_tilde2{3}(:),dual_tilde2{5}(:)];
        epi_h{2} = vec - sigma2*Proj_epih2(sigma2,vec);
        
        % Update variables
        dual_var2{1} = reshape(epi_h{1}(:,1),Nx,Ny);
        dual_var2{4} = reshape(epi_h{1}(:,2),Nx,Ny);
        
        dual_var2{2} = reshape(epi_h{2}(:,1),Nx,Ny);
        dual_var2{3} = reshape(epi_h{2}(:,2),Nx,Ny);
        dual_var2{5} = reshape(epi_h{2}(:,3),Nx,Ny);
        
        for i=1:P+2
            g2{i} = dual_var2{i};
        end
    end
    
    tm = toc(tm);

    % Check criteria
    criteria_check;
    
    %% Reweighting procedure
    
    check_cond = rm <= param_algo.rel_stop_crit && ((method == 3 && pol_thresh == 50) || method ~= 3)...
        && ((proj_l2 == 1 && min_res(t) <= (tol*epsilon{n_test})) || (proj_l2 ~= 1 && t>10));
    
    
    if proj_l2 && rw < param_algo.num_rw && (check_cond || t-t_prev >= param_algo.tmax_rw)
        
        rw = rw + 1;
        fprintf('weighted l1 regularised case, iter = %e\n',rw);
        disp('***********************')
        fprintf('Seed %i\n',n_test);
        fprintf('Iter %i\n',t);
        fprintf('L1 norm              = %e\n',l1_norm(t))
        fprintf('Residual             = %e\n',residual(t))
        fprintf('Obj function value   = %e\n',obj(t))
        fprintf('Rel sol change = %e\n', rm)
        fprintf('Residual = %e\n',min_res(t))
        fprintf('Polarization constraint not satisfied:%d \n',c_f(t));
        fprintf('Polarization constraint not satisfied after threshold:%d \n',p_thresh(t));
        fprintf('SNR_I = %e, SNR_Q = %e, SNR_U = %e \n',snr_v{1}(t), snr_v{2}(t), snr_v{3}(t));
        
        disp('***********************')
        
        for i =1:P
            comp_weight;
            eta{i} = eta_o(i)*weights{i};
        end
        
        % Weighted L1 problem
        St_im_sol{rw} = St_im; % Solution for each reweighting iteration
        snr_v_sol{rw} = snr_v;
        t_sol{rw} = t;
        mse_sol{rw} = mse;
        nrmse_sol{rw} = nrmse;
        p_thresh_sol{rw} = p_thresh;
        cf_sol{rw} = c_f;
        
        if rw == 1 && proj_l2 == 1
            method = param_algo.method_change;
            param_algo.rel_stop_crit = param_algo.rel_stop_crit_change;
            tol = param_algo.tol_rw;
            parameters_algo;
        end
        
    elseif rw == param_algo.num_rw && check_cond
        break
    elseif ~strcmp(dict,'SARA') && check_cond
            break 
    end
    
end

disp('***********************')
fprintf('Seed %i\n',n_test);
fprintf('Iter %i\n',t);
fprintf('L1 norm              = %e\n',l1_norm(t))
fprintf('Residual             = %e\n',residual(t))
fprintf('Obj function value   = %e\n',obj(t))
fprintf('Rel sol change = %e\n', rm);
fprintf('Residual = %e\n',min_res(t));
fprintf('Relative obj change = %e\n', rel_obj(t));
fprintf('Polarization constraint not satisfied:%d \n',c_f(t));
fprintf('Polarization constraint not satisfied after threshold:%d \n',p_thresh(t));
fprintf('SNR_I = %e, SNR_Q = %e, SNR_U = %e \n',snr_v{1}(t), snr_v{2}(t), snr_v{3}(t));

disp('***********************')



clear c_f
clear obj
clear residual
clear l1_norm
clear min_res

    
            
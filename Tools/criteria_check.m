% Check criteria
rm = 0;
l1_norm(t) = 0;
residual(t) = 0;
min_res(t) = 0;

for i = 1:P
    if method ~= 1
        %     l1_norm(t) = l1_norm(t) + eta(i)*sum(abs(Psitw{i}(St_im{i})));
        temp = Psitw{i}(St_im{i});
        temp1 = temp;
        if strcmp(dict,'SARA')
            l1_norm(t) = l1_norm(t) + sum(eta{i}.*abs(temp));
        end
        if strcmp(dict,'TV')
            u = temp(1:Nx*Ny);
            v = temp(Nx*Ny+1:end);
            temp = sqrt(u.^2 + v.^2);
            l1_norm(t) = l1_norm(t) + sum(eta(i).*abs(temp));
        end
        if strcmp(dict,'Dirac+TV')
            temp = temp1+ temp;
            l1_norm(t) = l1_norm(t) + sum(eta(i).*abs(temp));
        end
        
    end
    
    rel_sol{i}(t) = norm(St_im{i}(:) - St_im_old{i}(:))/norm(St_im_old{i}(:));
    rm = max(rel_sol{i}(t),rm);
    snr_v{i}(t) = calculate_snr(im{i},St_im{i});
    [mse{i}(t),nrmse{i}(t)] = calculate_mse(im{i},St_im{i});
    
    
    %% l2 Projection checking
    
    if bright ~= 1
        res_im{i} = Gw*A(St_im{i});
        residual(t) = residual(t) + 0.5*sum((abs(res_im{i}-y{1,n_test}{i,1})).^2);
        res_norm(i) = norm(res_im{i} - y{1,n_test}{i,1},2);
        min_res(t) = max(min_res(t),res_norm(i) - epsilon{n_test});
    end
end

if bright == 1
    
    s_m = [St_im{1}(:),St_im{2}(:), St_im{3}(:)];
    bb = s_m*L;
    p_s = zeros(length(yb{n_test,1}{1,1}),4);
    for j = 1:4
        Bb{j} = reshape(bb(:,j), Ny, Nx);
        p_s(:,j) = Gw*A(Bb{j}); % To simulate the operation of phi(S)
        
        if sep == 1
            res_norm(j) = norm(p_s(:,j) - yb{n_test,1}{j,1},2);
            min_res(t) = max(min_res(t), res_norm(j) - epsilon{n_test});
            residual(t) = residual(t) + 0.5*res_norm(j);
        end
    end
    
    if sep~=1
        res_norm = norm(p_s(:) - y_conc(:), 2);
        min_res(t) = res_norm - epsilon{n_test};
        residual(t) = 0.5*res_norm;
    end
    
    clear s_m;
    clear bb;
    clear p_s;
    clear res_norm;
end

if method ~= 1
    if proj_l2 == 1
        obj(t) = l1_norm(t);
    else
        obj(t) = l1_norm(t) + residual(t) ;
    end
else
    obj(t) = residual(t) ;
end



%**************************************************************************
% Checking constraints
%**************************************************************************
if P>1
    P_abs = sqrt(abs(St_im{2}).^2 + abs(St_im{3}).^2);
    [c_f(t),~] = size(find(P_abs - St_im{1} > 0));
    c_m(t) = max(max(P_abs-St_im{1}));
    c_mn(t) = max(max((P_abs-St_im{1})./St_im{1}));
    
    if method == 3
        
        % To check the constraint for epi_h1
        [c_h1(t),~] = size(find(-St_im{1}>St_im{4}));
        
        % To check the constraint for epi_h2
        [c_h2(t),~] = size(find(P_abs-St_im{5} > 1e-5));
       
        if (param.verbose == 3)
            if c_h1(t) > 0
                fprintf('Epi_h1 not satisfied:%d \n',c_h1(t));
            end
            if c_h2(t) > 0
                fprintf('Epi_h2 not satisfied:%d \n',c_h2(t));
            end
            if c_f(t) > 0
                fprintf('Polarization constraint not satisfied:%d \n',c_f(t));
            end
            
            figure(102)
            subplot 221, plot(c_h1), title('Count of unsatisfied epi h1')
            subplot 222, plot(c_h2), title('Count of unsatisfied epi h2')
            subplot 223, plot(c_f), title('Count of unsatisfied total polarization constraint')
            subplot 224, plot(c_m), title('Max of norm(Q,U)-I')
        end
    end
end




%% Stopping criterion and logs

% log
if (param.verbose == 1) && rem(t,1000) == 0
    fprintf('Iter %i\n',t);
    fprintf('L1 norm              = %e\n',l1_norm(t))
    fprintf('Residual             = %e\n',residual(t))
    fprintf('Obj function value   = %e\n',obj(t))
    %         fprintf('Rel sol norm change  = %e\n',rel_sol)
end

if (param.verbose == 2 && rem(t,1000) == 0)
    figure(100)
    subplot 331, semilogy(obj), title('Obj function')
    subplot 334, plot(snr_v{1}), title('SNR- I')
    subplot 335, plot(snr_v{2}), title('SNR- Q')
    subplot 336, plot(snr_v{3}), title('SNR- U')
    subplot 337, semilogy(rel_sol{1}), title('rel norm- I')
    subplot 338, semilogy(rel_sol{2}), title('rel norm- Q')
    subplot 339, semilogy(rel_sol{3}), title('rel norm- U')
    
    
    figure(101)
    subplot 231, imagesc(im{1}), title('True I'), axis image, axis off
    subplot 232, imagesc(im{2}), title('True Q'), axis image, axis off
    subplot 233, imagesc(im{3}), title('True U'), axis image, axis off
    subplot 234, imagesc(St_im{1}), title('Reconstructed I'), axis image, axis off
    subplot 235, imagesc(St_im{2}), title('Reconstructed Q'), axis image, axis off
    subplot 236, imagesc(St_im{3}), title('Reconstructed U'), axis image, axis off
  
    pause(0.1)
end

count_pol_thresh;
p_thresh(t) = pol_thresh;

    
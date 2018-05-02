% Set convergence parameters for the algorithm

sigma1 = param_algo.Eta1/(2);
sigma2 = param_algo.Eta2/2;
sigma3 = param_algo.Eta3/(2*evl);

verbosity = 0;

if strcmp(param.dict,'Dirac+TV')
    param_algo.Eta4 = 5*10^3;
    tv_op_norm = op_norm(Psit_tv,Psi_tv,[Ny, Nx], 1e-4, 200, verbosity);
    sigma4 = param_algo.Eta4/(2*tv_op_norm);
    nu = [4e-2,2e-3,5e-3];  % Regularization parameter for TV norm
    %  nu = [4e-1,2e-2,5e-2];  % Regularization parameter for TV norm
end


if proj_l2 == 1
    switch(method)
        case 1
            tau = 2/param_algo.Eta3; % step-size for gradient descent
        case 2
            tau = 2/(param_algo.Eta1+param_algo.Eta3);
            if strcmp(dict,'Dirac+TV')
                tau = 1.99/(param_algo.Eta1+param_algo.Eta3+param_algo.Eta4);
            end
        case 3
            tau = 2/(param_algo.Eta1+param_algo.Eta2+param_algo.Eta3);
            if strcmp(dict,'Dirac+TV')
                tau = 2/(evl+param_algo.Eta1+param_algo.Eta2+param_algo.Eta3+param_algo.Eta4);
            end
    end
    
else
    switch(method)
        case 1
            tau = 1.99/evl; % step-size for gradient descent
        case 2
            tau = 1.99/(evl+param_algo.Eta1);
            if strcmp(dict,'Dirac+TV')
                tau = 1.99/(evl+param_algo.Eta1+param_algo.Eta4);
            end
        case 3
            tau = 1.99/(evl+param_algo.Eta1+param_algo.Eta2);
            if strcmp(dict,'Dirac+TV')
                tau = 1.99/(evl+param_algo.Eta1+param_algo.Eta2+param_algo.Eta4);
            end
    end
end


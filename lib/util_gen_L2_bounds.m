function [epsilont, epsilonts, epsilon, epsilons] = util_gen_L2_bounds(y0, input_snr, s, equal_bounds)
% generates the input data

if ~exist('s', 'var')
    s = [2 2;
         0 0;
         0 0
         1 1];
else
    if size(s, 2) == 1
        s = ones(size(s, 1), 2)*s;
    end
    if size(s, 1) == 1
        s = [s; 0 0; 0 0; 1 1;];
    end
end
if ~exist('equal_bounds', 'var')
    equal_bounds = 0;
end


R = length(y0);
Nm = numel(cell2mat(y0));
normy0 = norm(cell2mat(y0));

% add Gaussian i.i.d. noise
sigma_noise = 10^(-input_snr/20) * normy0/sqrt(Nm);


if s(4, 1) == 1
    s1 = s(1, 1);
    % estimate L2 ball parameter
    epsilon = sqrt(Nm + s1*sqrt(2*Nm)) * sigma_noise;
    epsilont = cell(R,1);
    for q = 1:R
        % this produces a global bound which is greater than the mean by
        % ~sqrt(mean(length(y{:})) (if equal length)
        epsilont{q} = sqrt(size(y0{q}, 1) + s1*sqrt(2*size(y0{q}, 1))) * sigma_noise;
    end
end

if s(4, 2) == 1
    s2 = s(1, 2);
    % estimate L2 ball parameter
    epsilons = sqrt(Nm + s2*sqrt(2*Nm)) * sigma_noise;
    epsilonts = cell(R,1);
    for q = 1:R
        % this produces a global bound which is greater than the mean by
        % ~sqrt(mean(length(y{:})) (if equal length)
        epsilonts{q} = sqrt(size(y0{q}, 1) + s2*sqrt(2*size(y0{q}, 1))) * sigma_noise;
    end
end

if s(4, 1) == 2
    p1 = 1-s(2, 1);
    % estimate L2 ball parameter
    epsilon = sqrt(chi2inv(p1, Nm)) * sigma_noise;
    epsilont = cell(R,1);
    for q = 1:R
        % this produces a global bound which is greater
        epsilont{q} = sqrt(chi2inv(p1, size(y0{q}, 1))) * sigma_noise;
    end
end

if s(4, 2) == 2
    p2 = 1-s(2, 2);
    % estimate L2 ball parameter
    epsilons = sqrt(chi2inv(p2, Nm)) * sigma_noise;
    epsilonts = cell(R,1);
    for q = 1:R
       % this produces a global bound which is greater
        epsilonts{q} = sqrt(chi2inv(p2, size(y0{q}, 1))) * sigma_noise;
    end
end

if s(4, 2) == 3
    % estimate L2 ball parameter
    epsilons = epsilon * s(3, 2);
    epsilonts = cell(R,1);
    for q = 1:R
       % this produces a global bound which is greater
        epsilonts{q} = epsilont{q} * s(3, 2);
    end
end


if equal_bounds
    n1 = norm(cell2mat(epsilont))/epsilon;
    n2 = norm(cell2mat(epsilonts))/epsilons;
    for q = 1:R
        epsilont{q} = epsilont{q} / n1;
    end

    for q = 1:R
        epsilonts{q} = epsilonts{q} / n2;
    end
end

end


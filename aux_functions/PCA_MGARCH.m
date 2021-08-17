function [Hest, s2up_est] = PCA_MGARCH(dataset, K, N)
% INPUTS:
%   dataset      - A N by p matrix of zero mean data
%   K            - Number of common factors
%   C
% OUTPUTS:
%   Hest:       One-step-ahead of the conditional covariance matrix
%   Hone:       One-step-ahead of the conditional covariance matrix of the
%               common shocks
opt.disp = 0; 
univariteOptions = optimset('fminunc');
univariteOptions.Display = 'none';
univariteOptions.LargeScale = 'off';

[V, MM] = eigs(cov(dataset), K,'LM',opt);
uu = dataset*V;
R = V;
idioest = dataset-uu*R';
[~, H_one] = DCC_full_normal(uu);

for i = 1:N
    [pari, ~, hidio(:,i)] = tarch(idioest(:,i),1,0,1,[],[],[],univariteOptions);
    h_one(i) = pari(1) + pari(2)*idioest(end,i)^2 + pari(3)*hidio(end,i);
end
Hest = R*H_one*R'+diag(hidio(end,:)); 
end


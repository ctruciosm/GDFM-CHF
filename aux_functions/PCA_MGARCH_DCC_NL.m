function [Hest, s2up_est] = PCA_MGARCH_DCC_NL(dataset, K, N)
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

H_idio = DCCcov_fore(idioest,0);

Hest = R*H_one*R'+H_idio; 
end


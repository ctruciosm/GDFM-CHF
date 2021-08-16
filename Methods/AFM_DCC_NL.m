function [Hest, unc_Hu] = AFM_DCC_NL(dataset, K, N)
% Translation from zPOET.R (Mincheva et al. POET R package)
% INPUTS:
%   dataset      - A T (sample size) by N (assets) matrix of zero mean data
%   K            - Number of common factors
%   N            - Number of assets (time series)
% OUTPUTS:
%   SigmaU:     estimated p by p covariance matrix of y_t
%   SigmaY:     estimated p by p covariance matrix of u_t


opt.disp = 0; 
[V, MM] = eigs(cov(dataset), K,'LM',opt);
uu = dataset*V;
R = V;
idioest = dataset-uu*R';
unc_Hu = cov(uu);  
H_idio = DCCcov_fore(idioest,0);
Hest = R*unc_Hu*R'+H_idio; 
end


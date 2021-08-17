%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Estima??o de modelos GARCH multivariados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
warning off

addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/QuEST_v027'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/DCC_NL06'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/mfe-toolbox-master'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/Data'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/Methods'))


data = importdata('retornos_APP3_matlab.txt');

[T N] = size(data);
W = 750;
WR = T - 750;  

sigma_full = zeros(WR,N*N);

univariteOptions = optimset('fminunc');
univariteOptions.Display = 'none';
univariteOptions.LargeScale = 'off';

for l = 1:WR
    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));
    display(['Estimation OGARCH ',num2str(l),' out of ',num2str(WR)])
    [par,Ht,weights,pc,htMat] = o_mvgarch2(datatemp,N,1,0,1,[],univariteOptions) ;
    omega_est = par(1:3:N*3-2); 
    alpha_est = par(2:3:N*3-1);
    beta_est = par(3:3:N*3);

    s2ep_est = omega_est' + alpha_est'.*pc(end,:).^2 + beta_est'.*htMat(end,:);
    Haux = weights'*diag(s2ep_est)*weights;
    Hone = 0.5*(Haux+Haux');
    sigma_full(l,:) = Hone(:);
end

save('H_OGARCH.txt', 'sigma_full', '-ASCII');

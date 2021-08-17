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
L = T;
WR = L - 750;   

sigma_full = zeros(WR,N*N);

for l = 1:WR
    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));
    display(['Estimation DCC ',num2str(l),' out of ',num2str(WR)])
    [par, ~, H_u,Qt,intercept,stdData]= dcc_normal(datatemp,[],1,0,1,1,0,1,2,'3-stage','Diagonal');
    Q_one = intercept + par(end-1)*stdData(:,:,end) + par(end)*Qt(:,:,end);
    q_one = sqrt(diag(Q_one));
    R_one = Q_one./ (q_one*q_one');
    omega_est = par(1:3:N*3-2); 
    alpha_est = par(2:3:N*3-1);
    beta_est = par(3:3:N*3);
    D_one = zeros(N,N);
    for jj = 1:N
        D_one(jj,jj) = sqrt(omega_est(jj) + alpha_est(jj)*datatemp(end,jj)^2 + beta_est(jj)*H_u(jj,jj,end));
    end
    H_one = D_one*R_one*D_one;
    Hone = 0.5*(H_one+H_one');
    sigma_full(l,:) = Hone(:);
end
save('sigma_DCC.txt', 'sigma_full', '-ASCII');


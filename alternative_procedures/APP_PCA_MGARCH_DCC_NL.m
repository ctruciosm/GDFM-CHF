%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Unconditional Covariance matrix with non-linear shrinkage
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

sigma_NL = zeros(WR,N*N);


for l = 1:WR
    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));
    display(['Estimation PCA_MGARCH_DCC_NL APP3 ',num2str(l),' out of ',num2str(WR)])

    AUXHt = PCA_MGARCH_DCC_NL(datatemp, 8, N);
    Ht = 0.5*(AUXHt+AUXHt');
    sigma_NL(l,:) = Ht(:);
end
          
save('H_PCA_MGARCH_DCC_NL.txt', 'sigma_NL', '-ASCII');







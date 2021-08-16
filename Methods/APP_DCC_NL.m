%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Unconditional Covariance matrix with non-linear shrinkage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
warning off

addpath(genpath('/home/alunos/10/ra109078/QuEST_v027'))
addpath(genpath('/home/alunos/10/ra109078/DCC_NL06'))
addpath(genpath('/home/alunos/10/ra109078/GDFM_MC/MFE'))

%addpath(genpath('/Volumes/CTRUCIOS_SD/GDFM-CHF/NL/QuEST_v027'))
%addpath(genpath('/Volumes/CTRUCIOS_SD/GDFM-CHF/DCC-NL/DCC_NL06'))
%addpath(genpath('/Volumes/CTRUCIOS_SD/VaR_GDFM_CHF/CovarianceEstimationCodes/GDFM_VaR/MFE'))


data = importdata('retornos_APP3_matlab.txt');

[T N] = size(data);
W = 750;
L = T;
WR = L - 750;

sigma_NL = zeros(WR,N*N);


for l = 1:WR
    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));
    display(['Estimation DCC_NL APP3 ',num2str(l),' out of ',num2str(WR)])

    AUXHt = DCCcov_fore(datatemp,0);
    Ht = 0.5*(AUXHt+AUXHt');
    sigma_NL(l,:) = Ht(:);
end
          
save('H_DCC_NL.txt', 'sigma_NL', '-ASCII');







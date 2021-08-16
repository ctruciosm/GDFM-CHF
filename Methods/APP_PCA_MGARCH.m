%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Unconditional Covariance matrix with non-linear shrinkage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
warning off

addpath(genpath('/home/alunos/10/ra109078/QuEST_v027'))
addpath(genpath('/home/alunos/10/ra109078/DCC_NL06'))
addpath(genpath('/home/alunos/10/ra109078/GDFM_VaR/MFE'))
addpath(genpath('/home/alunos/10/ra109078/PCA_MGARCH'))
addpath(genpath('/home/alunos/10/ra109078/GDFM_VaR/Estimadores'))



% Macbook path
%addpath(genpath('/Volumes/CTRUCIOS_SD/Research/GDFM-CHF/NL/QuEST_v027'))
%addpath(genpath('/Volumes/CTRUCIOS_SD/Research/GDFM-CHF/DCC-NL/DCC_NL06'))
%addpath(genpath('/Volumes/CTRUCIOS_SD/Research/VaR_GDFM_CHF/CovarianceEstimationCodes/GDFM_VaR/MFE'))
%addpath(genpath('/Volumes/CTRUCIOS_SD/Research/GDFM-CHF/PCA_MGARCH'))
%addpath(genpath('/Volumes/CTRUCIOS_SD/Research/VaR_GDFM_CHF/CovarianceEstimationCodes/GDFM_VaR/Estimadores'))

% Windows path 
%addpath(genpath('D:\Carlos\Research\GDFM-CHF\QuEST_v027'))
%addpath(genpath('D:\Carlos\Research\GDFM-CHF\DCC_NL06'))
%addpath(genpath('D:\Carlos\Research\GDFM-CHF\GDFM_VaR\MFE'))
%addpath(genpath('D:\Carlos\Research\GDFM-CHF\PCA_MGARCH'))
%addpath(genpath('D:\Carlos\Research\GDFM-CHF\GDFM_VaR\Estimadores'))



data = importdata('retornos_APP3_matlab.txt');

[T N] = size(data);
W = 750;
L = T;
WR = L - 750;

sigma_NL = zeros(WR,N*N);


for l = 1:WR
    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));
    display(['Estimation PCA_MGARCH APP3 ',num2str(l),' out of ',num2str(WR)])

    AUXHt = PCA_MGARCH(datatemp, 8, N);
    Ht = 0.5*(AUXHt+AUXHt');
    sigma_NL(l,:) = Ht(:);
end
          
save('H_PCA_MGARCH.txt', 'sigma_NL', '-ASCII');







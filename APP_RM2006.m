%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Estima??????o de modelos GARCH multivariados
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
    display(['Estimation RM2006 APP3',num2str(l),' out of ',num2str(WR)])

    [Htfull, paramsDCC]=riskmetrics2006_fore(datatemp);
    Haux = Htfull(:,:,end);
    Hone = 0.5*(Haux+Haux');
    sigma_full(l,:) = Hone(:);
end

save('H_RM2006.txt', 'sigma_full', '-ASCII');







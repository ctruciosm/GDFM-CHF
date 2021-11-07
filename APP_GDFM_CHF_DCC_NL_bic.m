%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% UGDFM-CHF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
warning off

addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/QuEST_v027'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/DCC_NL06'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/mfe-toolbox-master'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/aux_functions'))
addpath(genpath('/Users/ctruciosm/Desktop/GDFM-CHF/Data'))


data = importdata('retornos_APP3_matlab.txt');

[T N] = size(data);
W = 750;
L = T;
WR = L - 750;

sigma_full = zeros(WR,N*N);

l = 1;
AUX_data = data(l:l+750-1,:);
datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));

[T N] = size(datatemp);
q_max = 10; 
nmin = round(3*N/4);
q_aux = HL2(datatemp,q_max,2,nmin,'p1')
q = q_aux(2);    
% end numbers of common shocks
K = 10;                 % Lag B(L)
k = 5;                  % max Lags for the VAR
m = floor(sqrt(T));     % Lag Spectral density matrix
nrepli = 30;            % number os permutations


for l = 1:WR
    randn('state',l);
    rand('state',l);
    unifrnd('state',l);

    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));

    display(['Estimation GDFM_CHF_DCC_NL APP3 ',num2str(l),' out of ',num2str(WR)])

    
    [chi, CL, v] = fhlz_nstd_p_bic(datatemp(1:end,:),q+1,k,m,K,1:q,nrepli);
    idioest = datatemp(k+1:end,:)-chi;
        
    [~, H_one] = DCC_full_normal(v);  
    Hidio = DCCcov_fore(idioest,0);
    
    Haux = CL(:,:,1)*H_one*CL(:,:,1)'+Hidio;
    Hone = 0.5*(Haux+Haux');
    sigma_full(l,:) = Hone(:);
end
          
save('H_GDFM_CHF_DCC_NL_normal_bic.txt', 'sigma_full', '-ASCII');




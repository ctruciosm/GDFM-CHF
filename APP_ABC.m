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

% Number of common shocks and common factors in the first panel
AUX_data = data(1:1+W-1,:);
mu = mean(AUX_data);
datatemp = bsxfun(@minus,AUX_data,mu);

nmin = round(3*N/4);
q_aux = HL2(datatemp,q_max,2,nmin,'p1')
q = q_aux(2);
nfactors = BN_crit2(datatemp,8);
% end numbers of common shocks and common factors

   

K = 10;                 % Lag B(L)
m = floor(sqrt(T));     % Lag Spectral density matrix
k = 1; 

univariteOptions = optimset('fminunc');
univariteOptions.Display = 'none';
univariteOptions.LargeScale = 'off';
offset = 0;


for l = 1:WR
    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));
    display(['Estimation ABC APP3 ',num2str(l),' out of ',num2str(WR)])
    [CL, chi, uest,lambda,G]  = ML_DfmRawImp_noWW(datatemp, q, nfactors, k, K);
    idioest = datatemp-chi;
    [par, ~, H_u,Qt,intercept,stdData]= dcc_normal(uest,[],1,0,1,1,0,1,2,'3-stage','none');
    Q_one = intercept + par(end-1)*stdData(:,:,end) + par(end)*Qt(:,:,end);
    q_one = sqrt(diag(Q_one));
    R_one = Q_one./ (q_one*q_one');
    omega_est = [par(1) par(4) par(7)];
    alpha_est = [par(2) par(5) par(8)];
    beta_est = [par(3) par(6) par(9)];
    D_one = zeros(q,q);
    for jj = 1:q
        D_one(jj,jj) = omega_est(jj) + alpha_est(jj)*uest(end,jj)^2 + beta_est(jj)*H_u(jj,jj,end);
    end
    H_one = D_one*R_one*D_one;

    for i =1:N
        [paridio(:,i), ~, hidio(:,i)] = tarch(idioest(:,i),1,0,1,'NORMAL',[],[],univariteOptions);
        s2ep_est(i) = paridio(1,i) + paridio(2,i)*idioest(end,i)^2 + paridio(3,i)*hidio(end,i);
    end

    Haux = lambda*G*H_one*G'*lambda'+diag(s2ep_est);
    Hone = 0.5*(Haux+Haux');
    sigma_full(l,:) = Hone(:);
               
end
               
save('sigma_FHZG.txt', 'sigma_full', '-ASCII');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%% Unconditional Covariance matrix with non-linear shrinkage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
warning off

addpath(genpath('/home/alunos/10/ra109078/QuEST_v027'))
addpath(genpath('/home/alunos/10/ra109078/DCC_NL06'))
addpath(genpath('/home/alunos/10/ra109078/GDFM_MC/MFE'))
addpath(genpath('/home/alunos/10/ra109078/GDFM_VaR'))


data = importdata('retornos_APP3_matlab.txt');

[T N] = size(data);
W = 750;
L = T;
WR = L - 750;

sigma_full = zeros(WR,N*N);

% Number of common shocks in the first panel
AUX_data = data(1:1+W-1,:);
mu = mean(AUX_data);
datatemp = bsxfun(@minus,AUX_data,mu);

q_max = 10;
nmin = round(3*N/4);
q_aux = HL2(datatemp,q_max,2,nmin,'p1')
q = q_aux(2);
% end numbers of common shocks

K = 10;                 % Lag B(L)
m = floor(sqrt(T));     % Lag Spectral density matrix
k = 1;                  % Lags for the VAR
nrepli = 30;            % number os permutations

for l = 1012:WR
    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));
    display(['Estimation GDFM_CHF_DCC_NL_BEKK APP3 ',num2str(l),' out of ',num2str(WR)])
    
    [chi, CL, v] = fhlz_nstd_p(datatemp(1:end,:),q+1,k,m,K,1:q,nrepli);
    idioest = datatemp(k+1:end,:)-chi;

    [parameters,~,Ht] = bekk(v,[],1,0,1,'Full');
    [C,A,G,B] = bekk_parameter_transform(parameters,1,0,1,q,'Full');
    temp = v(end,:)'*v(end,:);
    H_one = C + A'*temp*A + B'*Ht(:,:,end)*B;

    Hidio = DCCcov_fore(idioest,0);
    
    Haux = CL(:,:,1)*H_one*CL(:,:,1)'+Hidio;
    Hone = 0.5*(Haux+Haux');
    sigma_full(l,:) = Hone(:);
end
          
save('H_GDFM_CHF_DCC_NL_BEKK.txt', 'sigma_full', '-ASCII');







%%%%%%%%%%%%%%%%%%%
%%%%  2GDFM4V  %%%%
%%%%%%%%%%%%%%%%%%%

clear all
clc
warning off

addpath(genpath('/home/alunos/10/ra109078/BH_replication_gdfm_vol'))
addpath(genpath('/home/alunos/10/ra109078/QuEST_v027'))
addpath(genpath('/home/alunos/10/ra109078/DCC_NL06'))
addpath(genpath('/home/alunos/10/ra109078/GDFM_MC/MFE'))
addpath(genpath('/home/alunos/10/ra109078/GDFM_VaR'))

data = importdata('retornos_APP3_matlab.txt');

[T N] = size(data);
W = 750;
L = T;
WR = L - 750;
sigma_NL = zeros(WR,N*N);

%%%%
opt.disp = 0; 
univariteOptions = optimset('fminunc');
univariteOptions.Display = 'none';
univariteOptions.LargeScale = 'off';

% Begin numbers of common shocks (level and volatilities)
%   computed using the first window panel 
AUX_data = data(1:1+W-1,:);
mu = mean(AUX_data);
datatemp = bsxfun(@minus,AUX_data,mu);
q_max = 10;
nmin = round(3*N/4);
q_aux = HL2(datatemp,q_max,2,nmin,'p1')
q = q_aux(2);
qq = q+1
k = 1;
idvar = [1:q];
nlagsimp = 10;  
nrepli = 30;  
w = floor(sqrt(T));
[x,CL,u,C1,ee] = fhlz_nstd(datatemp,q+1,k,w,nlagsimp,idvar,nrepli);
z=datatemp(k+1:end,:)-x;
ee = ee(k+1:end,:);          
for j=1:N
        [A,vv(:,j)]=ML_VAR(z(:,j),k,0);   
end
ss = log(ee.^2);  
oo = log(vv.^2);  
q_s_aux = HL2(ss,q_max,2,nmin,'p1')
q_s = q_s_aux(2);
q_o_aux = HL2(oo,q_max,2,nmin,'p1')
q_o = q_o_aux(2);
% end numbers of common shocks (level and volatilities)
Ht = zeros(N,N);

for l = 1:WR
    AUX_data = data(l:l+750-1,:);
    datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));
    display(['Estimation 2GDFM4V APP3 ',num2str(l),' out of ',num2str(WR)])
    
    % First GDFM
    [x,CL,u,C1,ee] = fhlz_nstd(datatemp,qq,k,w,nlagsimp,idvar,nrepli);
    z=datatemp(k+1:end,:)-x;
    ee = ee(k+1:end,:);                                                     % Eq (2.3.5) in BH17        
    for j=1:N
        [A,vv(:,j)]=ML_VAR(z(:,j),k,0);                                     % Eq. (2.5) in BH17
    end
    ee2 = (ee.^2);
    vv2 = (vv.^2);
    ss = log(ee2);                                                          % Eq (2.7) in BH17
    oo = log(vv2);                                                          % Eq (2.7) in BH17
    M_ss = mean(ss);
    M_oo = mean(oo);
    

    [chi_ss_fhlz,~,u_ss,~,~] = fhlz_nstd(ss,q_s+1,1,w,nlagsimp,[1:q_s],nrepli);  % Eq (2.8) BH17
    [chi_oo_fhlz,~,u_oo,~,~] = fhlz_nstd(oo,q_o+1,1,w,nlagsimp,[1:q_o],nrepli);  % Eq (2.9) BH17        
    
    varepsilon_ss_fhlz = ss(k+1:end,:)-chi_ss_fhlz;
    varepsilon_oo_fhlz = oo(k+1:end,:)-chi_oo_fhlz;
    
    
    % Apply GARCH(1,1) in each one
    exp_chi_ss_fhlz = exp(chi_ss_fhlz);
    exp_varepsilon_ss_fhlz = exp(varepsilon_ss_fhlz);
    exp_chi_oo_fhlz = exp(chi_oo_fhlz);
    exp_varepsilon_oo_fhlz = exp(varepsilon_oo_fhlz);
    
    % We will have the one-step-ahead conditional variance
    % Calculate Eq. (3.16), (3.17) and finally (3.18) which is my
    % one-step-ahead cond. covariance for each series 
    
    
    for ii = 1:N
        [par_exp_chi_ss_fhlz, ~, w_exp_chi_ss_fhlz(:,ii)] = tarch(exp_chi_ss_fhlz(:,ii),1,0,1,[],[],[],univariteOptions);
        [par_exp_varepsilon_ss_fhlz, ~, h_exp_varepsilon_ss_fhlz(:,ii)] = tarch(exp_varepsilon_ss_fhlz(:,ii),1,0,1,[],[],[],univariteOptions);
        [par_exp_chi_oo_fhlz, ~, w_exp_chi_oo_fhlz(:,ii)] = tarch(exp_chi_oo_fhlz(:,ii),1,0,1,[],[],[],univariteOptions);
        [par_exp_varepsilon_oo_fhlz, ~, h_exp_varepsilon_oo_fhlz(:,ii)] = tarch(exp_varepsilon_oo_fhlz(:,ii),1,0,1,[],[],[],univariteOptions);
    
         w_one_exp_chi_ss_fhlz(ii) = par_exp_chi_ss_fhlz(1) + par_exp_chi_ss_fhlz(2)*exp_chi_ss_fhlz(end,ii)^2 + par_exp_chi_ss_fhlz(3)*w_exp_chi_ss_fhlz(end,ii);
         h_one_exp_varepsilon_ss_fhlz(ii) = par_exp_varepsilon_ss_fhlz(1) + par_exp_varepsilon_ss_fhlz(2)*exp_varepsilon_ss_fhlz(end,ii)^2 + par_exp_varepsilon_ss_fhlz(3)*h_exp_varepsilon_ss_fhlz(end,ii);
         w_one_exp_chi_oo_fhlz(ii) = par_exp_chi_oo_fhlz(1) + par_exp_chi_oo_fhlz(2)*exp_chi_oo_fhlz(end,ii)^2 + par_exp_chi_oo_fhlz(3)*w_exp_chi_oo_fhlz(end,ii);
         h_one_exp_varepsilon_oo_fhlz(ii) = par_exp_varepsilon_oo_fhlz(1) + par_exp_varepsilon_oo_fhlz(2)*exp_varepsilon_oo_fhlz(end,ii)^2 + par_exp_varepsilon_oo_fhlz(3)*h_exp_varepsilon_oo_fhlz(end,ii);
    
         Ht(ii,ii) = w_one_exp_chi_ss_fhlz(ii)*h_one_exp_varepsilon_ss_fhlz(ii)*exp(M_ss(1,ii)) + w_one_exp_chi_oo_fhlz(ii)*h_one_exp_varepsilon_oo_fhlz(ii)*exp(M_oo(1,ii));
    end
    sigma_NL(l,:) = Ht(:);
end
          
save('H_2GDFM4V.txt', 'sigma_NL', '-ASCII');




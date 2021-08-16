function [SigmaU, SigmaY, factors, loadings] = POET(Y, K, C, thres, matrix)
% Translation from zPOET.R (Mincheva et al. POET R package)
% INPUTS:
%   Y            - A p by N matrix of zero mean data
%   K            - Number of common factors
%   C
%   thres = 'soft'
%   matrix = 'cor'
% OUTPUTS:
%   SigmaU:     estimated p by p covariance matrix of y_t
%   SigmaY:     estimated p by p covariance matrix of u_t
%   factors:    estimated unobservable factors in a K by T matrix form
%   loadings:   estimated factor loadings in a p by K matrix form


addpath(genpath('/Volumes/CTRUCIOS_SD/GDFM-CHF/NL/QuEST_v027'))
addpath(genpath('/Volumes/CTRUCIOS_SD/GDFM-CHF/DCC-NL/DCC_NL06'))
addpath(genpath('/Volumes/CTRUCIOS_SD/VaR_GDFM_CHF/CovarianceEstimationCodes/GDFM_VaR/MFE'))
K = 8;

data = importdata('retornos_APP3_matlab.txt');
AUX_data = data(1:750,:);
datatemp = bsxfun(@minus,AUX_data,mean(AUX_data));

Y = transpose(datatemp);


matrix = 'cor';
thres = 'soft';
C = 0.5;

[p n] = size(Y);

if K>0
    [V,Dd] = eig(transpose(Y)*Y);
    [W,I] = sort(Dd(:));
    F = sqrt(n) * V(:,1:K);
    LamPCA = Y * F/n;
    uhat = Y - LamPCA * transpose(F);
    Lowrank = LamPCA * transpose(LamPCA);
    rate = 1/sqrt(p) + sqrt((log(p))/n);
else
    uhat = Y;
    rate = sqrt((log(p))/n);
    Lowrank = zeros(p);
end

SuPCA = uhat * transpose(uhat)/n;
SuDiag = diag(diag(SuPCA));


if matrix =='cor'
    R = inv(SuDiag^(1/2)) * SuPCA * inv(SuDiag^(1/2));
end
if matrix  =='vad'
    R=SuPCA;
end

uu = zeros(p,p,n);
roottheta = zeros(p,p);
lambda = zeros(p,p);

for i = 1:p                  %adaptive threshold
    for j = 1:p
        uu(i,j,:) = uhat(i,:).*uhat(j,:);
        roottheta(i,j) = std(uu(i,j,:));
        lambda(i,j) = roottheta(i,j)*rate*C;
        lambda(j,i) = lambda(i,j);
    end
end

Rthresh = zeros(p,p);
if thres=='soft'
    for i = 1:p
        for j = 1:i
            if abs(R(i,j)) < lambda(i,j) && j < i
                Rthresh(i,j) = 0;
            else
                if j == i
                    Rthresh(i,j) = R(i,j);
                else
                    Rthresh(i,j) = sign(R(i,j))*(abs(R(i,j))-lambda(i,j));
                end
            end
            Rthresh(j,i) = Rthresh(i,j);
        end
    end
end
            
SigmaU = zeros(p,p);

if matrix =='cor'
    SigmaU = SuDiag^(1/2)*Rthresh*SuDiag^(1/2);
end
if matrix =='vad'
    SigmaU = Rthresh;
end

SigmaY = SigmaU + Lowrank;
factors = transpose(F);
loadings = LamPCA;

end








    

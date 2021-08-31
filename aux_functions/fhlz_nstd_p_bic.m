% FHLZ - **** **** *****
% This code is a fork of the one in the Matteo Barigozzi website
% Some minor modifications were made
% X  - data matrix
% qq - should correspond to (q+1) in the paper
% k  - number of minimum lag for the polynomial A(L)
% w  - number of lag for the covariance matrix estimation
% nlagsimp
%
%   H(L)*x(t)=R*v(t)+H(L)*xsi(t)
%   chi(t)=C(L)v(t) --> H(L)*chi(t)=R*v(t) 
%                                   if C(L) is zeroless
%           Aj(L)*chi(t)=Rj*vj(t) for j=1:m
%
%
%X= datatemp(1:end,:);
%qq = q+1;
%w = m;
%nlagsimp = K;
%idvar=1:q;
%nrepli=30;
function [cstar,CL,v] = fhlz_nstd_p_bic(X,qq,k,w,nlagsimp,idvar,nrepli)

T = size(X,1);
q = length(idvar);                                                         

sigma = std(X);
mu = mean(X);
z = (X - ones(T,1)*mu);

covmat = ML_cestimate(z,q,w);                                              
[n , o , h] = size(covmat);
H =(h+1)/2;
m = floor(n/qq);                                                            

for h = 1:nrepli
    imp = nan*ones(n,q,nlagsimp);
    riord = randperm(n);                                                    
    while not(isempty(intersect(idvar ,riord(qq*m+1:end))))
        riord = randperm(n);                                                
    end    
    x = z(:,riord);                                                         
    covmatrix = covmat(riord,riord,:);                                         
    [kkbic, kkaic, kkhq] = irfbic(covmatrix,q,k,T);
    [HL GL GL1 y]=StaticRepresentationbic(x,covmatrix,qq,m,nlagsimp,T,kkaic,H);     
    S = diag(std(y)); yy = ML_center(y)*(S^-1);  
    opt.disp = 0; 
    [V, MM] = eigs(cov(yy), q,'LM',opt);                                   
    M = diag(sqrt(diag(MM)));                                                                                                  % v(t)=sqrt(inv(Lambda))P'w(t)
    uu = yy*V*inv(M);  
    for lag = 1 : nlagsimp; 
        BBB(:, :, lag) = GL(:, :, lag)*S*V*M;          
    end;
    imp(riord(1:qq*m),:,:) = BBB;                                           
    tau = size(uu,1);
    [beta(:,:,:,h), u(:,:,h)] = DfmCholIdent(imp,idvar,uu(tau-T+k+1:end,:));
       
end

CL = nanmean(beta,4);
v = nanmean(u,3);

vv=[zeros(nlagsimp-1,q); v];
cstar=zeros(n,T-k);
for ii=1:T-k;
    for jj=1:nlagsimp;
        cstar(:,ii)=cstar(:,ii)+CL(:,:,jj)*vv(ii+nlagsimp-jj,:)';
    end;
end;
cstar=cstar';
end


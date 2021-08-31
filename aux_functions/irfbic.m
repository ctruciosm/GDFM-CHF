function [bic, aic, hq] = irfbic(covmatrix, q, K, T)
[n, ~, ~] = size(covmatrix);
qq = q  + 1;
s = floor(n/qq);
bic = zeros(s,K);
aic = zeros(s,K);
hq  = zeros(s,K);
 for k = 1:K
     for ss = 1:s-1
% compute relevant variances and covariances
[A, B] = yule(covmatrix((ss-1)*qq+1:ss*qq, (ss-1)*qq+1:ss*qq, :), k-1); 
% compute the VAR coefficients
C = A\B;
% compute the variance of the residuals
Ss = A(1:qq, 1:qq) - B'*C;
% compute the bic 
bic(ss,k) = log(det(Ss)) + (qq^2*k)*log(T)/T;
aic(ss,k) = log(det(Ss)) + (1/T)*2*k*qq^2;
hq(ss,k) = log(det(Ss)) + (qq^2*k+qq)*2*log(log(T))/T;
     end
a = qq*(s-1);
nlast = n-a;
% repeat the steps above for the last group of variables which does not have size $q+1$
[A, B] = yule(covmatrix(a+1:n, a+1:n, :), k-1);  
C = A\B;
Ss = A(1:nlast, 1:nlast) - B'*C;
bic(s,k) = log(det(Ss)) + (nlast^2*k)*log(T)/T;
hq(s,k) = log(det(Ss)) + (nlast^2*k+nlast)*2*log(log(T))/T;
 end
% minimize bic and cbic 
 [~, bic] = min(bic,[],2);
 [~, aic] = min(aic,[],2);
 [~, hq] = min(hq,[],2);
end

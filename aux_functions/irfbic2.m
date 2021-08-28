% rispetto ad irfbic qui nei casi in cui l'OLS e' mal
% condizionato l'output e' bic=0, aic=0, hq=0
% n.b: covmatrix is a 3-dimensional object having the lags on the third
% dimension (the central lag is lag 0); q is the no. of shocks, K is the
% maximum VAR lag, T is the no. of time obs. pval is the probability of the
% right tail
function [bic, aic, hq] = irfbic2(covmatrix, q, K, T)
[n, ~, ~] = size(covmatrix);
qq  = q  + 1;
s   = floor(n/qq);
bic = zeros(s,K);
aic = zeros(s,K);
hq  = zeros(s,K);
badcond=0;
for k = 1:K
    for ss = 1:s-1
        % compute relevant variances and covariances
        [A, B] = yule(covmatrix((ss-1)*qq+1:ss*qq, (ss-1)*qq+1:ss*qq, :), k-1);
        tolA=1e-9;
        if rcond(A)>tolA
            % compute the VAR coefficients
            C = A\B;
            % compute the variance of the residuals
            Ss = A(1:qq, 1:qq) - B'*C;
            % compute the bic
            bic(ss,k) = log(det(Ss)) + (qq^2*k)*log(T)/T;
            aic(ss,k) = log(det(Ss)) + (1/T)*2*k*qq^2;
            hq(ss,k) = log(det(Ss)) + (qq^2*k+qq)*2*log(log(T))/T;
        else
            badcond = 1;
        end
        
    end
    if ~badcond
        a = qq*(s-1);
        nlast = n-a;
        % repeat the steps above for the last group of variables which does not have size $q+1$
        [A, B] = yule(covmatrix(a+1:n, a+1:n, :), k-1);
        tolA=1e-9;
        if rcond(A)>tolA
            C = A\B;
            Ss = A(1:nlast, 1:nlast) - B'*C;
            bic(s,k) = log(det(Ss)) + (nlast^2*k)*log(T)/T;
            aic(s,k) = log(det(Ss)) + (1/T)*2*k*nlast^2;
            hq(s,k) = log(det(Ss)) + (nlast^2*k+nlast)*2*log(log(T))/T;
        else
            badcond = 1;
        end
    end
end

if ~badcond
    % minimize bic and cbic
    [~, bic] = min(bic,[],2);
    [~, aic] = min(aic,[],2);
    [~, hq]  = min(hq,[],2);
else
    [bic, aic, hq] = deal(0,0,0);
end

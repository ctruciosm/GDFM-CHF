function [ll,lls,Rt] = dcc_likelihood2(parameters,data,m,n,R,N,backCast,stage,composite,isJoint,isInference,gScale,univariate)
% Likelihood for estimation of scalar DCC(m,n) multivarate volatility models 
% with with GARCH(p,q) conditional variances
%
% USAGE:
%  [LL,LLS,RT,TB] = dcc_likelihood2(PARAMETERS,DATA,DATAASYM,M,L,N,R,N,BACKCAST,BACKCASTASYM,STAGE,COMPOSITE,ISJOINT,ISINFERENCE,GSCALE,UNIVARIATE)
%
% INPUTS:
%   PARAMETERS   - Vector of ADCC parameters, and possibly volatility and intercept parameters, 
%                    depending on other inputs
%   DATA         - A T by K matrix of zero mean residuals -OR-
%                    K by K by T array of covariance estimators (e.g. realized covariance)
%   M            - Order of symmetric innovations in DCC model
%   N            - Order of lagged correlation in DCC model
%   R            - K by K correlation matrix of standardized data
%   N            - K by K matrix mean of asymmetric standardized data outer products
%   BACKCAST     - K by K  matrix to use for back casting symetric terms
%   STAGE        - Integer, either 2 or 3 indicating whether 2-stage ro 3-stage estimator is being used
%   COMPOSITE    - Integer, one of 0 (None, use QMLE), 1 (Use diagonal composite) or 2 (full composite)
%   ISJOINT      - Boolean indicating whether PARAMETERS includes volatility parameters
%   ISINFERENCE  - Boolean indicating whether likelihood is used for making inference, in which case
%                    no transformations are made to parameters.
%   GSCALE       - K by 1 vector used in 2-stage to scale the intercept.  See DCC.
%   UNIVARIATE   - Cell array of structures containing information needed to compute volatilities.
%
% OUTPUTS:
%   LL           - The log likelihood evaluated at the PARAMETERS
%   LLS          - A T by 1 vector of log-likelihoods
%   RT           - A [K K T] dimension matrix of conditional correlations
%
% COMMENTS:
%
% See also DCC

% Copyright: Kevin Sheppard
% kevin.sheppard@economics.ox.ac.uk
% Revision: 1    Date: 4/13/2012
% Modified by Olivier Ledoit to save some memory
% assumes no asymmetric terms

% Modes:
% 1-stage: Univariate, R (corr vech), DCC (N-scale ignored)
% 2-stage (Estimation): R (transformed corr vech), DCC  (N-scale treated as constants)
% 2-stage (Inference): Univariate, R (corr vech), DCC (N-scale ignored)
% 3-stage (Estimation): DCC
% 3-stage (Inferecne): Univariate, R (corr vech), N (vech), DCC

if ndims(data)==2
    [T,k] = size(data);
    data2d = data;
    stdData = zeros(k,k,T);
    for t=1:T
        stdData(:,:,t) = data2d(t,:)'*data2d(t,:);
    end
elseif ndims(data)==3
    [k,~,T]=size(data);
    stdData=data;
    clear data
else
    error('DATA must be either a K by T matrix of a K by K by T 3-dimensional array.')
end

offset = 0;
% Parse Parameters
if stage==1 || isJoint
    count = 0;
    for i=1:k
        u = univariate{i};
        count = count + u.p+u.o+u.q+1;
    end
    garchParameters = parameters(1:count);
    offset = offset + count;
    computeVol = true;
else
    computeVol = false;
end
if stage<=2 || isJoint
    count = k*(k-1)/2;
    R = parameters(offset + (1:count));
    offset = offset + count;
end
a = parameters(offset + (1:m));
b = parameters(offset + (m+1:m+n));
if isempty(b)
    b=0;
end

% Compute volatilities
H = ones(T,k);
if computeVol
    H = dcc_reconstruct_variance(garchParameters,univariate);
    for t=1:T
        h = sqrt(H(t,:));
        stdData(:,:,t) = stdData(:,:,t)./(h'*h);
    end
    logdetH = sum(log(H),2);
else
    logdetH = zeros(T,1);
end

% Transform R & N, if needed
if stage<=2 || isJoint
    % Transform R
    if isInference
        R = corr_ivech(R);
    else
        R = z2r(R);
    end
end

% Compute intercept
if stage==3
    intercept = R*max(1-sum(a)-sum(b),0);
else
    scale = (1-sum(a)-sum(b)) - gScale*sum(g);
    scale = sqrt(scale);
    intercept = R.*(scale*scale');
end

% Indices or constant, as needed
if composite == 0
    likconst = k*log(2*pi);
elseif composite == 1
    indices = [(1:k-1)' (2:k)'];
elseif composite == 2
    [i,j] = meshgrid(1:k);
    indices = [i(~triu(true(k))) j(~triu(true(k)))];
end

I = eye(k);
Qt = zeros(k,k,n+1);
Rt = zeros(k,k,T);
lls = zeros(T,1);
for t=1:T
    Qt(:,:,n+1) = intercept;
    for i = 1:m
        if (t-i)>0
            Qt(:,:,n+1) = Qt(:,:,n+1) + a(i)*stdData(:,:,t-i);
        else
            Qt(:,:,n+1) = Qt(:,:,n+1) + a(i)*backCast;
        end
    end
    for i = 1:n
        if (t-i)>0
            Qt(:,:,n+1) = Qt(:,:,n+1) + b(i)*Qt(:,:,n+1-i);
        else
            Qt(:,:,n+1) = Qt(:,:,n+1) + b(i)*backCast;
        end
    end
    q = sqrt(diag(Qt(:,:,n+1)));
    Rt(:,:,t) = Qt(:,:,n+1)./ (q*q');
    if composite == 0
        lls(t) = 0.5*(likconst + logdetH(t) + log(det(Rt(:,:,t))) + sum(diag((Rt(:,:,t)\I)*stdData(:,:,t))));
    elseif composite
        S = (sqrt(H(t,:))'*sqrt(H(t,:))) .* Rt(:,:,t);
        lls(t) = composite_likelihood(S,stdData(:,:,t),indices);
    end
    Qt(:,:,1:n)=Qt(:,:,2:n+1);
end
ll = sum(lls);

if isnan(ll) || ~isreal(ll) || isinf(ll)
    ll = 1e7;
end

% compute number of bytes utilized in memory by Matlab in this function
if 0&(isequal(computer,'PCWIN')|isequal(computer,'PCWIN64'))
   userview=memory;
   disp(userview.MemUsedMATLAB);
end


function [sigmahat, Ht, parameters, ll]=DCCcov(x,demean)

% function [sigmahat, Ht, parameters, ll]=DCCcov(x,demean)
%
% Estimation of scalar DCC(1,1) multivarate volatility model with GARCH(1,1)
% conditional variances
%
% INPUTS:
%   x            - A T by K matrix of zero mean residuals
%   demean       - set to 1 if the data need to be demeaned first, 
%                  and 0 otherwise (optional, default is 1)
%
% OUTPUTS:
%   sigmahat     - K*K conditional covariance matrix matrix on last day of
%                  sample [it is equal to Ht(:,:,end)]
%   Ht           - A [K K T] dimension matrix of conditional covariances
%   PARAMETERS   - Estimated parameters.
%                    3-stage: [VOL(1) ... VOL(K) corr_vech(R)' vech(N)' alpha gamma beta]
%                    where VOL(j) is a 3-dimensional vector containing the parameters from
%                    volatility model j.
%   LL           - The log likelihood at the optimum
%
% COMMENTS:
%   The dynamics of the correlations in a 3-stage DCC model are:
%     Q(t) = R*(1-a-b) + a*e(t-1)'*e(t-1) + b*Q(t-1)
%
% same as dccstream11.m but with only straight DCC-NL model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    demean=1;
end
startingVals = [];
[T,k] = size(x);
% demean data if desired
if demean
    meanx=mean(x);
    x=x-meanx(ones(T,1),:);
end
stdData = zeros(k,k,T);
for t=1:T
   stdData(:,:,t) = x(t,:)'*x(t,:);
end
options = optimset('fmincon');
options.Display = 'none';
options.Diagnostics = 'off';
options.Algorithm = 'interior-point';
options.TolX=1e-6;
options.TolFun=1e-4;
optimset(options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Univariate volatility models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[H,univariate] = dcc_fit_variance(x,ones(k,1),zeros(k,1), ...
   ones(k,1),2*ones(k,1),[]);
for t=1:T
   h = sqrt(H(t,:));
   hh = h'*h;
   stdData(:,:,t) = stdData(:,:,t)./hh;
end

%%%%%%%%%%%%
% Back casts
%%%%%%%%%%%%
w = .06*.94.^(0:sqrt(T));
w = w'/sum(w);
backCast = zeros(k);
for i=1:length(w)
   backCast = backCast + w(i)*stdData(:,:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of the unconditional correlation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C = mean(stdData,3);
correl = sqrt(diag(C));
C = C ./ (correl*correl');
C = (C + C') ./ 2;

% use nonlinear shrinkage to estimate unconditional covariance matrix of
% devolatilized residuals
R = QuESTimate(x./sqrt(H),0,[1e-5 40]);
end
r = sqrt(diag(R));
R = R ./ (r*r');
R = (R + R') ./ 2;

% now estimate the two dynamic parameters of the DCC model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [.01 .03 .05 .1];
theta = [.99 .97 .95];
[a,theta] = ndgrid(a,theta);
parameters = unique([a(:) theta(:)-a(:)],'rows');
minLL= inf;
for i=1:size(parameters,1)
   ll = dcc_likestream(parameters(i,:),stdData,R,backCast);
   if ll<minLL
      startingVals = parameters(i,:);
      minLL = ll;
   end
end
a = startingVals(1);
b = startingVals(2);
startingVals = [a b];
LB = zeros(length(startingVals),1);
UB = ones(length(startingVals),1);
A = ones(1,length(startingVals));
b = .99998;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of the ARCH/GARCH parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parameters = fmincon(@dcc_likestream,startingVals,A,b,[],[],LB,UB,[], ...
   options,stdData,R,backCast);
clear stdData

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Covariances and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
garchParameters = [];
for i=1:k
   garchParameters  = [garchParameters univariate{i}.parameters'];
end
parameters = [garchParameters corr_vech(R)' parameters];

% compute log-likelihood at the optimum
[ll,~,Ht] = dcc_likelihood2(parameters,x,1,1,[],[],backCast,3,1,true,true,[],univariate);
ll = -ll;

% recombine conditional correlations with conditional variances
for t=1:T
   h = sqrt(H(t,:));
   Ht(:,:,t) = Ht(:,:,t).*(h'*h);
end
sigmahat=Ht(:,:,T);

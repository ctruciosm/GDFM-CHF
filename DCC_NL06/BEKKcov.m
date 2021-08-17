function [sigmahat, Ht, parameters, ll, intercept]=BEKKcov(x,demean)

% function  [sigmahat, Ht, parameters, ll, intercept]=BEKKcov(x,demean)
% 
% Estimation of covariance targeting scalar BEKK(1,0,1) multivarate
% volatility model using Contiguous Pairs Composite Likelihoods
%
% INPUTS:
%   X            - A T by K matrix of zero mean residuals
%   demean       - set to 1 if the data need to be demeaned first, 
%                  and 0 otherwise (optional, default is 1)
%
% OUTPUTS:
%   sigmahat     - K*K conditional covariance matrix matrix on last day of
%                  sample [it is equal to Ht(:,:,end)]
%   HT           - A [K K T] dimension matrix of conditional covariances
%   PARAMETERS   - Estimated parameters.  
%                   alpha beta
%   LL           - The log likelihood at the optimum
%   INTERCEPT    - K by K matrix containing the intercept computed from the unconditional variance
%                    and parameters   
%
% COMMENTS:
%    The conditional variance, H(t), of a scalar variance-targeting vech is modeled
%    as follows:
%
%      H(t) = (1-alpha-beta)*C + alpha*r_{t-1}'*r_{t-1} + beta*H(t-1)
%
% same as Zhao Zhao's vtstream.m function but only doing nonlinear
% shrinkage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Argument Checking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<2
    demean=1;
end

%data should be TxK, T>K
[t,k]=size(x);
% demean data if desired
if demean
    meanx=mean(x);
    x=x-meanx(ones(t,1),:);
end
data = zeros(k,k,t);
for i=1:t
    data(:,:,i) = x(i,:)'*x(i,:);
end
options = optimset('fmincon');
options.Display = 'none';
options.Diagnostics = 'off';
options.Algorithm = 'interior-point';
options.TolX=1e-6;
options.TolFun=1e-4;
optimset(options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the backCast
backCast = zeros(k);
tau = max(ceil(sqrt(t)),k);
weights = .06 * .94.^(0:tau);
weights = weights / sum(weights);
for i=1:tau
    backCast = backCast + weights(i) * data(:,:,i);
end

C0 = mean(data,3);
% change the effective sample size for shrinkage
tEff = 150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimation of the unconditional covariance matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replace C with nonlinear shrinkage C
C = QuESTimate(x,0,[1e-5 40]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting Values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = [.01 .03 .05 .1];
theta = [.99 .97 .95];
[a,theta] = ndgrid(a,theta);
parameters = unique([a(:) theta(:)-a(:)],'rows');
minLL= inf;
for i=1:size(parameters,1)
   %   ll = dcc_likelihood(parameters(i,:),stdData,[],1,0,1,R,[],backCast,[],3,1,false,false);
   ll = vt_likestream(parameters(i,:),data,C,backCast);
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
parameters = fmincon(@vt_likestream,startingVals,A,b,[],[],LB,UB,[], ...
   options,data,C,backCast);

[ll,~,Ht]=scalar_vt_vech_likelihood(parameters,data,[],1,0,1,C,[],[],backCast,[],false,1,false);
ll=-ll;
sigmahat=Ht(:,:,t);
clear data

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the intercept
%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = parameters(1);
beta = parameters(2);
intercept = C*(1-alpha-beta);

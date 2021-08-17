% BEKK_NLdemo.m - demonstrate how to estimate the scalar BEKK-NL model
% originally based on DCC_NLdemo.m but adapted for the BEKK-NL model
% details of the BEKK-NL model are in Section 7 of the article 
% entitled "Large Dynamic Covariance Matrices"  in the Journal of 
% Business and Economics Statistics, authored by Robert F. Engle, 
% Olivier Ledoit, and Michael Wolf

clear
clear global

% add path to QuEST toolbox [would be different on your computer]
QuESTpath='D:\Research\QuEST\QuEST_v027\';
addpath(QuESTpath)

% add path to direct kernel estimator [would be different on your computer]
DIRECTpath='D:\Research\Direct\Public_Code\';
addpath(DIRECTpath)

% add path to Oxford MFE toolbox [would be different on your computer]
mfepath='D:\Research\kevinsheppard-mfe_toolbox-45462370b5f0\';
CURRENT_DIRECTORY=pwd;
cd(mfepath)
addToPath('-silent')
cd(CURRENT_DIRECTORY)
clear CURRENT_DIRECTORY

% load test data set (if too big for your RAM, cut off some columns)
load DCC_NLdata.mat x
x=x(:,1:100);
x=double(x);
[t,n]=size(x);

% estimate BEKK-NL model
disp(['Starting to estimate BEKK-NL model at ' datestr(now)])
sigmahat=BEKKcov(x);
disp(['Finished estimating BEKK-NL model at ' datestr(now)])

% display a sample of the results for verification purposes
disp(sigmahat(end-2:end,end-2:end)*100)

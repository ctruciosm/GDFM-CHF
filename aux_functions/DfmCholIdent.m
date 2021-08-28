function [imp, sh] = DfmCholIdent(rawimp, variables,rsh)
% rawimp = imp;
% variables = idvar;
% rsh= uu(tau-T+K+1:end,:)


sh=[];
q = length(variables);
B0 = rawimp(variables,1:q,1);
C = chol(B0*B0')';
H = inv(B0)*C;
k = size(rawimp,3);
for j =  1:k
    imp(:,:,j) = rawimp(:,:,j)*H;
end
if nargin == 3
    sh = rsh*H;
end

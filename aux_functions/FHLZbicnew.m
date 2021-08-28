
function [BBB, u, AA, R, fore,mindeta]= FHLZbicnew(y, covmatrix, q, k, nlagsimp)

sd = std(y);
x = y*diag(sd.^(-1));
[ T, n] = size(x);
qq = q  + 1;
s = floor(n/qq);
if isscalar(k)
  k = k*ones(s,1);
end
kk = max(k);
% inizialize matrices
AA = zeros(n,n,kk+1);
BB = zeros(n,n,nlagsimp);
BBB = zeros(n,q,nlagsimp);
for ss = 1:s-1
% compute relevant variances and covariances
[A, B] = yule(covmatrix((ss-1)*qq+1:ss*qq, (ss-1)*qq+1:ss*qq, :), k(ss)-1); 
% compute the determinant
detA(ss) = det(A);
% compute coefficients for each small VAR
C = A\B;
CC = zeros(qq,qq,k(ss) + 1);
CC(:,:,1) = eye(qq);
CC(:,:,2:end) = - reshape(C', qq, qq, k(ss));
% compute the big VAR 
AA((ss-1)*qq+1:qq*ss,(ss-1)*qq+1:qq*ss, 1:k(ss)+1) = CC;
% compute piecewise the VAR inverse
BB((ss-1)*qq+1:qq*ss,(ss-1)*qq+1:qq*ss, :) = invertepolynomialmatrix(CC, nlagsimp);
end



a = qq*(s-1);
nlast = n-a;
% repeat the steps above for the last group of variables which does not have size $q+1$
[A, B] = yule(covmatrix(a+1:n, a+1:n, :), k(s)-1);  
C = A\B;
detA(s)= det(A);
mindeta = min(detA);
CC = zeros(nlast,nlast,k(s) + 1);
CC(:,:,1) = eye(nlast);
CC(:,:,2:k(s) + 1) = - reshape(C', nlast, nlast, k(s));
AA(a+1:n, a+1:n, 1:k(s)+1) = CC;
BB(a+1:n, a+1:n, :) = invertepolynomialmatrix(CC, nlagsimp);
% compute R
AAA = reshape(AA,n, n*(kk+1));
xx = zeros(T - kk, n*(kk+1));
for j = 1:kk+1
     xx(:,(j-1)*n+1:n*j) =  x(kk+2-j:T+1-j , :);
end
z = xx*AAA';
norm = 1./std(z);
z = z*diag(norm);
S = cov(z);
opt.disp = 0;
[V, MM] = eigs(S, q,'LM',opt);
 M = diag(sqrt(diag(MM)));
 R = V*M;
 R = diag(1./norm)*R;
u = z*V/M;

for lag = 1 : nlagsimp 
     BBB(:, :, lag) = diag(sd)*BB(:, :, lag)*R;
end
if nargout >= 5
   fore = prediction(BBB, u, 1);
end
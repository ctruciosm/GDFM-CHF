% inversion of a matrix of polynomials in the lag operator
function inverse = InvPolMatrix(poly,nlags); 
n = size(poly,1); k = size(poly,3) - 1;
for s = 1:k+1; newpoly(:,:,s)=inv(poly(:,:,1))*poly(:,:,s); end;            % Guarantees that poly(:,:,1)=eye(n)
polynomialmatrix = - newpoly(:,:,2:k+1);
A = zeros(n*k,n*k);
A(n+1:n*k,1:n*(k-1)) = eye(n*(k-1));
for j = 1:k
   A(1:n,(j-1)*n+1:j*n) = polynomialmatrix(:,:,j);
end
inverse = zeros(n,n,nlags);
D = eye(n*k);
for j = 1:nlags 
   inverse(:,:,j) = D(1:n,1:n)*inv(poly(:,:,1));
   D = A*D;
end
%%% =================================================================== %%%
end
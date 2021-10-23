
function  E = subr_1(x, M)
%eigv = subr_1(x(1:T,1:N), M); 
%x = center(x);
% define some useful quantities
[T,N] = size(x);
w = M;
W = 2*w+1;
% compute covariances
SS = zeros(W,N*N);
for k = 1:w+1,
     Sk = center(x(k:T,:))'*center(x(1:T+1-k,:))/(T-k);
     S_k = Sk';
     % the x(:) creates a column vector from the matrix x columnwise
     SS(w-k+2,:) = S_k(:)';  
     SS(w+k,:) = Sk(:)';
end


% compute the spectral matrix in w points (S)
% [0:2*pi/W:4*pi*w/W]
%Factor = exp(-sqrt(-1)*(-w:w)'*(0:2*pi/W:4*pi*w/W));


freq = [0:2*pi/W:pi];
Factor = exp(-sqrt(-1)*(-w:w)'*freq);

ww = 1 - abs([-w:w])/(w+1);

% multiply first line of SS by ww(1)
% multiply second line of SS by ww(2)

S = diag(ww)*SS(1:W,:); 

% inverse to the (:) operation
% S = reshape(S'*Factor,N,N*W);

S = reshape(S'*Factor,N,N*(w+1));

% compute the eigenvalues 
eigenvalues = zeros(N,w+1);
D = eig(S(:,1:N));
eigenvalues(:,1) = flipud(sort(real(D)));

for j = 1:w,
   D = eig(S(:,j*N+1:(j+1)*N));
   eigenvalues(:,1+j) = flipud(sort(real(D)));
end
%eigenvalues
%[eigenvalues(:,1)  eigenvalues(:,2:jj+1)*2]
E = [eigenvalues(:,1)  eigenvalues(:,2:w+1)*2]*ones(w+1,1)/(2*w+1);
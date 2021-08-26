function [S0, S1] = yule(covmatrix, k)
[n, ~, h] = size(covmatrix);
H = (h+1)/2;
covmat = covmatrix(:,:, H - k - 1:H + k + 1);
S0 = zeros(n*(k+1), n*(k+1));
striscia = reshape(covmat, n, n*(2*k+3));
for j = 1:k+1
    S0((j-1)*n+1:j*n,:) = striscia(:,(k+2-j)*n + 1:(2*k+3-j)*n);
end
  S1 = striscia(:,(k+2)*n + 1:(2*k+3)*n)'; 
end
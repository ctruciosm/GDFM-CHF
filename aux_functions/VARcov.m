%%% =================================================================== %%% 
% VARcov estimate a VAR from autocovariance of the variables
% C - autoregressive coefficients 
% B - Moving average Represetation
% u - residuals
function [C CC B u]=VARcov(x,covmatrix,T,k,qq,mm,H,nlagsimp);

A = zeros(qq*k,qq*k); B = zeros(qq*k,qq); xx = zeros(T-k,qq*k);
for j = 1:k
    for i = 1:k;    
        A((j-1)*qq+1:qq*j,(i-1)*qq+1:qq*i) = covmatrix((mm-1)*qq+1:mm*qq, (mm-1)*qq+1:mm*qq,H + i - j);   
    end    
    B((j-1)*qq+1:qq*j,:) = covmatrix((mm-1)*qq+1:mm*qq, (mm-1)*qq+1:mm*qq,H - j);     
    xx(:,(j-1)*qq+1:qq*j) =  x(k+1-j:T-j ,(mm-1)*qq+1:mm*qq);
end


C = inv(A)*B;                                                               % Aj(L)
u = x(k+1:T , (mm-1)*qq+1:mm*qq) - xx*C;                                    % Residual wj(t)
CC(:,:,1) = eye(qq); CC(:,:,2:k + 1) = - reshape(C',qq,qq,k);               % Reshaping Aj(L) s.t. it is invertible
B = InvPolMatrix(CC,nlagsimp);                                              % Building Gj(L)= inv(Aj(L))
end




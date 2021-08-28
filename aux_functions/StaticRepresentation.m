function [HL GL GL1 w]=StaticRepresentation(x,covmatrix,qq,m,nlagsimp,T,k,H);

GL = zeros(qq*m,qq*m,nlagsimp); HL = zeros(qq*m,qq*m,k+1); w = zeros(T-k,qq*m);

for mm = 1:m;                
    [C CC B u]=VARcov(x,covmatrix,T,k,qq,mm,H,nlagsimp);                    % Estimation of Aj(L) - Gj(L) - vj(t)    
    HL((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm,:)=CC;                           % Building the H(L) Matrix
    w(:,(mm-1)*qq+1:qq*mm) = u;                                             % Residual w(t)   
    GL((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm,:)=B;                            % Building G(L)=inv(H(L))
    GL1((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm,:)=inv(sum(CC,3));              % long run impact matrix   
end
end
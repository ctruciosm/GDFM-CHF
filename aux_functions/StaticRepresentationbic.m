%%% =================================================================== %%%
%%% =================================================================== %%%
function [HL GL GL1 w]=StaticRepresentationbic(x,covmatrix,qq,m,nlagsimp,T,k,H);

kk = max(k);
HL = zeros(qq*m,qq*m,kk+1);
GL = zeros(qq*m,qq*m,nlagsimp); 

for mm = 1:m;                
    [C CC B u]=VARcov(x,covmatrix,T,k(mm),qq,mm,H,nlagsimp);                % Estimation of Aj(L) - Gj(L) - vj(t)    
    HL((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm, 1:k(mm)+1)=CC;                   % Building the H(L) Matrix
    GL((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm,:)=B;                            % Building G(L)=inv(H(L))
    GL1((mm-1)*qq+1:qq*mm,(mm-1)*qq+1:qq*mm,:)=inv(sum(CC,3));              % long run impact matrix   
end
AAA = reshape(HL,qq*m, qq*m*(kk+1));
xx = zeros(T - kk, qq*m*(kk+1));
for j = 1:kk+1
     xx(:,(j-1)*qq*m+1:qq*m*j) =  x(kk+2-j:T+1-j ,1:qq*m);
end
w = xx*AAA';
end
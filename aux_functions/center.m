function XC = center(X)

[T n] = size(X);
XC = X - ones(T,1)*mean(X); 
XC = XC./kron(ones(T,1),std(XC));
end

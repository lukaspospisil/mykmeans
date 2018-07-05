function L = get_L( X,C,gamma )

[n,T] = size(X);

L = norm(X - C*gamma,'fro')^2/(T*n);

end


function N = Neumann2DCurv(G, m, n, b)
    % G is the curvilinear gradient and b is the Neumann coeff.

    Bm = sparse(m+2, m+1);
    Bm(1, 1) = -b;
    Bm(end, end) = b;
   
    Bn = sparse(n+2, n+1);
    Bn(1, 1) = -b;
    Bn(end, end) = b;
   
    Im = sparse(m + 2, m);
    In = sparse(n + 2, n);
   
    Im(2:(m+2)-1, :) = speye(m, m);
    In(2:(n+2)-1, :) = speye(n, n);
   
    Bm = kron(In, Bm);
    Bn = kron(Bn, Im);
   
    N = [Bm Bn]*G;
   
    % Then L + N + (N - R) should not be singular
end
function N = Neumann3DCurv(G, m, n, o, b)
    % G is the curvilinear gradient and b is the Neumann coeff.

    Bm = sparse(m+2, m+1);
    Bm(1, 1) = -b;
    Bm(end, end) = b;
   
    Bn = sparse(n+2, n+1);
    Bn(1, 1) = -b;
    Bn(end, end) = b;
   
    Bo = sparse(o+2, o+1);
    Bo(1, 1) = -b;
    Bo(end, end) = b;
   
    Im = sparse(m+2, m);
    In = sparse(n+2, n);
    Io = sparse(o+2, o);
   
    Im(2:(m+2)-1, :) = speye(m, m);
    In(2:(n+2)-1, :) = speye(n, n);
    Io(2:(o+2)-1, :) = speye(o, o);
   
    Bm = kron(kron(Io, In), Bm);
    Bn = kron(kron(Io, Bn), Im);
    Bo = kron(kron(Bo, In), Im);
   
    N = [Bm Bn Bo]*G;
   
    % Then L + Ncurv + (Rrect - Ncurv) should not be singular
end
function [du, dv] = d2Dc(m, n, dx, dy, forward_x, forward_y)
    du = d2_dx(m, n, dx, forward_x);
    dv = d2_dy(m, n, dy, forward_y);
end

function d = d2_dx(m, n, dx, forward)
    d = d_dx(m, dx, forward);
    I = speye(n);    
    d = kron(I', d);
end

function d = d2_dy(m, n, dy, forward)
    d = d_dx(n, dy, forward);
    I = speye(m);    
    d = kron(d, I');
end

function d = d_dx(m, dx, forward)
    if forward
        d = spdiags([-ones(m, 1) ones(m, 1)], [0 1], m, m);
        d(end, end-1) = -1;
        d(end, end) = 1;
    else
        d = spdiags([-ones(m, 1) ones(m, 1)], [-1 0], m, m);
        d(1, 1) = -1;
        d(1, 2) = 1;
    end
    d = d/dx;
end
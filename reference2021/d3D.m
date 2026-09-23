function [du, dv, dw] = d3D(m, n, o, dx, dy, dz, forward_x, forward_y, forward_z)
    du = d3_dx(m, n, o, dx, forward_x);
    dv = d3_dy(m, n, o, dy, forward_y);
    dw = d3_dz(m, n, o, dz, forward_z);
end

function d = d3_dx(m, n, o, dx, forward)
    d = d_dx(m, dx, forward);
    Io = speye(o);
    In = speye(n);
    d = kron(kron(Io', In'), d);
end

function d = d3_dy(m, n, o, dy, forward)
    d = d_dx(n, dy, forward);
    Io = speye(o);
    Im = speye(m);
    d = kron(kron(Io', d), Im');
end

function d = d3_dz(m, n, o, dz, forward)
    d = d_dx(o, dz, forward);
    In = speye(n);
    Im = speye(m);
    d = kron(kron(d, In'), Im');
end

function d = d_dx(m, dx, forward)
    if forward
        d = spdiags([-ones(m+1, 1) ones(m+1, 1)], [0 1], m+1, m+1);
        d(end, end-1) = -1;
        d(end, end) = 1;
    else
        d = spdiags([-ones(m+1, 1) ones(m+1, 1)], [-1 0], m+1, m+1);
        d(1, 1) = -1;
        d(1, 2) = 1;
    end
    d = d/dx;
end


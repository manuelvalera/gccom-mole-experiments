function I = GI13(M, m, n, o, type)
% Interpolates a logical-gradient component from the staggered face set it
% is computed on, to the face set it is needed on, for grad3DCurv.
%
% Parameters:
%          M : the operator block to be interpolated (Ge, Gn or Gc)
%    m, n, o : number of CELLS along xi, eta, zeta
%       type : which shift to perform (see table below)
%
% Face layouts, flattened xi-fastest then eta then zeta:
%     x-faces : (m+1) x  n    x  o     xi on nodes
%     y-faces :  m    x (n+1) x  o     eta on nodes
%     z-faces :  m    x  n    x (o+1)  zeta on nodes
%
% Each shift is a tensor product of three 1-D operators, one per axis:
%     type    source -> target      xi      eta     zeta
%     'Gn'    y-face -> x-face      c->n    n->c     I
%     'Gc'    z-face -> x-face      c->n     I      n->c
%     'Ge'    x-face -> y-face      n->c    c->n     I
%     'Gcy'   z-face -> y-face       I      c->n    n->c
%     'Gee'   x-face -> z-face      n->c     I      c->n
%     'Gnn'   y-face -> z-face       I      n->c    c->n

    switch type
        case 'Gn',  Ax = P(m); Ay = Q(n); Az = speye(o);
        case 'Gc',  Ax = P(m); Ay = speye(n); Az = Q(o);
        case 'Ge',  Ax = Q(m); Ay = P(n); Az = speye(o);
        case 'Gcy', Ax = speye(m); Ay = P(n); Az = Q(o);
        case 'Gee', Ax = Q(m); Ay = speye(n); Az = P(o);
        case 'Gnn', Ax = speye(m); Ay = Q(n); Az = P(o);
        otherwise,  error('GI13:BadType', 'unknown type "%s"', type);
    end

    I = kron(Az, kron(Ay, Ax)) * M;
end

function A = Q(N)
% node -> center, N+1 values on nodes -> N values on centers. Midpoint average.
    A = spdiags(0.5*ones(N, 2), [0 1], N, N+1);
end

function A = P(N)
% center -> node, N values on centers -> N+1 values on nodes.
% Midpoint average in the interior, linear extrapolation at the two ends.
    A = spdiags(0.5*ones(N+1, 2), [-1 0], N+1, N);
    A(1, 1) = 1.5;  A(1, 2) = -0.5;              % extrapolate to node 0
    A(N+1, N) = 1.5;  A(N+1, N-1) = -0.5;        % extrapolate to node N
end

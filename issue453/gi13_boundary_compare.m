1;
% gi13_boundary_compare.m -- which boundary row should GI13's centre->node
% operator P use?  Run from the ROOT of the mole-fork checkout:
%
%   octave-cli path\to\gi13_boundary_compare.m
%
% The Copilot review on PR #466 points out that GI2 (the 2-D analogue) uses a
% three-point boundary row [1, .5, -.5], while the PR uses the two-point
% linear extrapolation [1.5, -.5].  Both are exact for linear fields and both
% are second order; they differ in the error constant.  At the boundary node,
% for f = x^2 with centres at h/2, 3h/2, 5h/2:
%
%     [1.5, -.5]       error -0.75 h^2
%     [1, .5, -.5]     error -1.75 h^2      (GI2)
%
% So the PR's row is ~2.3x more accurate at the boundary node, and GI2's is
% the convention the 2-D code already uses.  This script measures what that
% difference does to grad3DCurv itself -- interior faces, boundary faces, and
% the convergence order -- so the choice can be made on data.

src = fullfile(pwd, 'src', 'octave');
if ~exist(fullfile(src, 'grad3DCurv.m'), 'file')
  error('run this from the root of the mole checkout (no src/octave/grad3DCurv.m here)');
end
addpath(genpath(src));
warning('off', 'all');

function write_gi13(d, variant)
  if strcmp(variant, 'ours')
    brow = ["    A(1, 1) = 1.5;  A(1, 2) = -0.5;\n" ...
            "    A(N+1, N) = 1.5;  A(N+1, N-1) = -0.5;\n"];
  else
    brow = ["    A(1, 1) = 1;  A(1, 2) = 0.5;  A(1, 3) = -0.5;\n" ...
            "    A(N+1, N) = 1;  A(N+1, N-1) = 0.5;  A(N+1, N-2) = -0.5;\n"];
  end
  txt = ["function I = GI13(M, m, n, o, type)\n" ...
         "    switch type\n" ...
         "        case 'Gn',  Ax = P(m); Ay = Q(n); Az = speye(o);\n" ...
         "        case 'Gc',  Ax = P(m); Ay = speye(n); Az = Q(o);\n" ...
         "        case 'Ge',  Ax = Q(m); Ay = P(n); Az = speye(o);\n" ...
         "        case 'Gcy', Ax = speye(m); Ay = P(n); Az = Q(o);\n" ...
         "        case 'Gee', Ax = Q(m); Ay = speye(n); Az = P(o);\n" ...
         "        case 'Gnn', Ax = speye(m); Ay = Q(n); Az = P(o);\n" ...
         "        otherwise,  error('GI13:BadType', 'unknown type');\n" ...
         "    end\n" ...
         "    I = kron(Az, kron(Ay, Ax)) * M;\n" ...
         "end\n" ...
         "function A = Q(N)\n" ...
         "    A = spdiags(0.5*ones(N, 2), [0 1], N, N+1);\n" ...
         "end\n" ...
         "function A = P(N)\n" ...
         "    A = spdiags(0.5*ones(N+1, 2), [-1 0], N+1, N);\n" ...
         brow ...
         "end\n"];
  fid = fopen(fullfile(d, 'GI13.m'), 'w'); fprintf(fid, '%s', txt); fclose(fid);
end

function [eall, eint, ebnd] = grad_err(N, amp, k)
  [a, b, c] = meshgrid(linspace(0,1,N), linspace(0,1,N), linspace(0,1,N));
  X = a + amp*sin(2*pi*b); Y = b + amp*sin(2*pi*a);
  Z = c + 0.5*amp*sin(2*pi*a).*sin(2*pi*b);
  s = linspace(0,1,N); s = [0, s(1:end-1)+0.5/(N-1), 1];
  [p, q, t] = meshgrid(s, s, s);
  x = p + amp*sin(2*pi*q); y = q + amp*sin(2*pi*p);
  z = t + 0.5*amp*sin(2*pi*p).*sin(2*pi*q);
  f = reshape(permute(x.^2 + y.^2 + z.^2, [2 1 3]), [], 1);
  T = grad3DCurv(k, X, Y, Z)*f;
  cx = linspace(0,1,N); cx = cx(1:end-1) + 0.5/(N-1);
  [aa, bb, ~] = meshgrid(linspace(0,1,N), cx, cx);
  xf = aa + amp*sin(2*pi*bb);
  nx = N*(N-1)*(N-1);
  d = reshape(T(1:nx) - reshape(permute(2*xf, [2 1 3]), [], 1), N, N-1, N-1);
  eall = sqrt(mean(d(:).^2));
  bd = d([1 end], :, :);     ebnd = sqrt(mean(bd(:).^2));   % where P's boundary row acts
  in = d(2:end-1, :, :);     eint = sqrt(mean(in(:).^2));
end

function r = lin_err(N, amp)
  [a, b, c] = meshgrid(linspace(0,1,N), linspace(0,1,N), linspace(0,1,N));
  X = a + amp*sin(2*pi*b); Y = b + amp*sin(2*pi*a);
  Z = c + 0.5*amp*sin(2*pi*a).*sin(2*pi*b);
  s = linspace(0,1,N); s = [0, s(1:end-1)+0.5/(N-1), 1];
  [p, q, t] = meshgrid(s, s, s);
  x = p + amp*sin(2*pi*q); y = q + amp*sin(2*pi*p);
  z = t + 0.5*amp*sin(2*pi*p).*sin(2*pi*q);
  f = reshape(permute(x + 2*y + 3*z, [2 1 3]), [], 1);
  T = grad3DCurv(2, X, Y, Z)*f;
  r = sqrt(mean((T(1:N*(N-1)*(N-1)) - 1).^2));
end

tmp = tempname(); mkdir(tmp);
names = {'ours', 'gi2'};
labels = {'[1.5, -.5]   (PR as submitted)', '[1, .5, -.5]  (GI2 convention)'};
Ns = [13 21 33 49];
R = struct();

for v = 1:2
  d = fullfile(tmp, names{v}); mkdir(d); write_gi13(d, names{v});
  addpath(d, '-begin'); rehash(); clear GI13;
  w = which('GI13');
  if isempty(strfind(w, d))
    error('variant %s not picked up; which(GI13) = %s', names{v}, w);
  end
  printf('\n=== P boundary row %s ===\n', labels{v});
  printf('  %4s  %-12s %-6s %-12s %-6s %-12s %-6s\n', 'n', 'all faces', 'p', ...
         'interior', 'p', 'xi-boundary', 'p');
  E = zeros(numel(Ns), 3);
  for i = 1:numel(Ns)
    [E(i,1), E(i,2), E(i,3)] = grad_err(Ns(i), 0.10, 2);
    if i == 1
      printf('  %4d  %-12.4e %-6s %-12.4e %-6s %-12.4e %-6s\n', Ns(i), ...
             E(i,1), '-', E(i,2), '-', E(i,3), '-');
    else
      rr = log((Ns(i)-1)/(Ns(i-1)-1));
      o = log(E(i-1,:)./E(i,:))/rr;
      printf('  %4d  %-12.4e %-6.2f %-12.4e %-6.2f %-12.4e %-6.2f\n', Ns(i), ...
             E(i,1), o(1), E(i,2), o(2), E(i,3), o(3));
    end
  end
  R.(names{v}) = E;
  printf('  linear field f=x+2y+3z, n=21: rms %.4e\n', lin_err(21, 0.10));
  rmpath(d); clear GI13;
end

printf('\n=== ratio GI2 / ours (above 1 means the PR row is more accurate) ===\n');
printf('  %4s  %-10s %-10s %-10s\n', 'n', 'all', 'interior', 'boundary');
for i = 1:numel(Ns)
  q = R.gi2(i,:)./R.ours(i,:);
  printf('  %4d  %-10.3f %-10.3f %-10.3f\n', Ns(i), q(1), q(2), q(3));
end
printf('\nBoth rows are second order.  If the ratios are close to 1 the choice is\n');
printf('pure convention and matching GI2 costs nothing.  If the boundary ratio is\n');
printf('well above 1, matching GI2 trades accuracy for consistency -- say so in\n');
printf('the reply either way.\n');

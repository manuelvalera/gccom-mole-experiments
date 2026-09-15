% test_grad3DCurv_convergence.m
%
% Convergence study for the 3D curvilinear gradient and for L = D*G.
% Run once with MOLE's GI13.m on the path, then again with the patched
% GI13.m ahead of it, and compare the two tables.
%
% Works in Octave and MATLAB.

MOLE_SRC = '../src/matlab_octave';     % <-- adjust
addpath(MOLE_SRC);
% addpath('/path/to/patched');         % <-- uncomment to test the patch
                                       %     (addpath prepends, so this must
                                       %      come AFTER the MOLE_SRC line)

k   = 2;
amp = 0.10;                            % grid distortion

fprintf('GI13 in use: %s\n\n', which('GI13'));
fprintf('3D curvilinear grid, distortion amp = %.2f, k = %d\n', amp, k);
fprintf('field f = x^2 + y^2 + z^2   (grad_x = 2x, laplacian = 6)\n\n');
fprintf('%5s %13s %7s %13s %7s\n', 'n', 'grad rms', 'order', 'L=D*G rms', 'order');

pg = NaN; pl = NaN; pn = NaN;
for N = [13 21 33 49]
  [a0, b0, c0] = meshgrid(linspace(0,1,N), linspace(0,1,N), linspace(0,1,N));
  X = a0 + amp*sin(2*pi*b0);
  Y = b0 + amp*sin(2*pi*a0);
  Z = c0 + 0.5*amp*sin(2*pi*a0).*sin(2*pi*b0);
  sc = linspace(0,1,N);  sc = [0, sc(1:end-1) + 0.5/(N-1), 1];
  [a1, b1, c1] = meshgrid(sc, sc, sc);
  xc = a1 + amp*sin(2*pi*b1);
  yc = b1 + amp*sin(2*pi*a1);
  zc = c1 + 0.5*amp*sin(2*pi*a1).*sin(2*pi*b1);
  f = reshape(permute(xc.^2 + yc.^2 + zc.^2, [2 1 3]), [], 1);

  % --- gradient, x-component, against the exact 2x on the x-faces
  G  = grad3DCurv(k, X, Y, Z);
  T  = G*f;
  s  = linspace(0, 1, N);
  cc = s(1:end-1) + 0.5/(N-1);
  [a, b, ~] = meshgrid(s, cc, cc);
  xf = a + amp*sin(2*pi*b);
  gx = T(1:N*(N-1)*(N-1));
  eg = sqrt(mean((gx - reshape(permute(2*xf, [2 1 3]), [], 1)).^2));

  % --- Laplacian, interior only (no BC operator applied here)
  L = div3DCurv(k, X, Y, Z)*G;
  v = reshape(L*f, N+1, N+1, N+1);
  d = v(3:end-2, 3:end-2, 3:end-2) - 6;
  el = sqrt(mean(d(:).^2));

  if isnan(pg)
    og = '   -'; ol = '   -';
  else
    r  = log((N-1)/(pn-1));
    og = sprintf('%5.2f', log(pg/eg)/r);
    ol = sprintf('%5.2f', log(pl/el)/r);
  end
  fprintf('%5d %13.4e %7s %13.4e %7s\n', N, eg, og, el, ol);
  pg = eg; pl = el; pn = N;
end

% --- regression: a Cartesian grid must stay exact ------------------------
fprintf('\nregression, Cartesian grid (must stay at roundoff):\n');
for N = [11 17 25]
  [X, Y, Z]    = meshgrid(linspace(0,1,N), linspace(0,1,N), linspace(0,1,N));
  sc = linspace(0,1,N);  sc = [0, sc(1:end-1) + 0.5/(N-1), 1];
  [xc, yc, zc] = meshgrid(sc, sc, sc);
  f = reshape(permute(xc.^2 + yc.^2 + zc.^2, [2 1 3]), [], 1);
  L = div3DCurv(k, X, Y, Z)*grad3DCurv(k, X, Y, Z);
  v = reshape(L*f, N+1, N+1, N+1);
  d = v(3:end-2, 3:end-2, 3:end-2) - 6;
  fprintf('  n=%2d  rms %.3e  max %.3e\n', N, sqrt(mean(d(:).^2)), max(abs(d(:))));
end

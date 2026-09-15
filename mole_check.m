% mole_check.m -- all four MOLE curvilinear-operator checks in one file.
%
% USAGE
%   1. Edit MOLE_SRC below to point at your MOLE clone's matlab_octave folder.
%   2. octave mole_check.m        (or type  mole_check  at the Octave prompt)
%
% Needs nothing but Octave. Writes its own scratch files to a temp folder.

MOLE_SRC = 'C:/Users/dap21/gccom-mole-experiments-2026/mole/src/matlab_octave';

% ---------------------------------------------------------------- scratch
scratch = fullfile(tempdir(), 'mole_check');
gdir    = fullfile(scratch, 'grids');
sdir    = fullfile(gdir, 'seamount');
pdir    = fullfile(scratch, 'patch');
for d = {scratch, gdir, sdir, pdir}
  if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

  function wf(p, txt)
    f = fopen(p, 'w'); fputs(f, txt); fclose(f);
  end

wf(fullfile(sdir, 'stretch.m'), [ ...
"function xp = stretch(s, L, beta)\n" ...
"  D = L/2;\n" ...
"  A = (1/(2*beta))*log((1+(exp(beta)-1)*D/L)/(1+(exp(-beta)-1)*D/L));\n" ...
"  xp = D*(1 + sinh(beta*(s-A))/sinh(beta*A));\n" ...
"end\n"]);

wf(fullfile(sdir, 'betaOf.m'), [ ...
"function bb = betaOf(s, which)\n" ...
"  global BT BB\n" ...
"  if isempty(BT), BT=2.0; end\n" ...
"  if isempty(BB), BB=5.5; end\n" ...
"  if strcmp(which,'bot'), bb = BB*(1 - 0.35*sin(pi*s));\n" ...
"  else,                   bb = BT*(1 + 0.55*sin(pi*s));\n" ...
"  end\n" ...
"end\n"]);

wf(fullfile(sdir, 'bottom.m'), [ ...
"function XY = bottom(s)\n" ...
"  Ls=1000; a=0.5; b=8; x0=-1.8; x1=1.8; xc=0.7;\n" ...
"  x = x0 + (x1-x0)*stretch(s,1,betaOf(s,'bot'));\n" ...
"  XY = [x*Ls, Ls*(-1 + a*exp(-b*(x-xc)^2))];\n" ...
"end\n"]);

wf(fullfile(sdir, 'top.m'), [ ...
"function XY = top(s)\n" ...
"  Ls=1000; x0=-1.8; x1=1.8;\n" ...
"  x = x0 + (x1-x0)*stretch(s,1,betaOf(s,'top'));\n" ...
"  XY = [x*Ls, 0];\n" ...
"end\n"]);

wf(fullfile(sdir, 'left.m'), [ ...
"function XY = left(s)\n" ...
"  Ls=1000; x0=-1.8; a=0.5; b=8; xc=0.7;\n" ...
"  d = Ls*(-1 + a*exp(-b*(x0-xc)^2));\n" ...
"  XY = [(x0 + 0.12*sin(pi*s))*Ls, d*(1-s)];\n" ...
"end\n"]);

wf(fullfile(sdir, 'right.m'), [ ...
"function XY = right(s)\n" ...
"  Ls=1000; x1=1.8; a=0.5; b=8; xc=0.7;\n" ...
"  d = Ls*(-1 + a*exp(-b*(x1-xc)^2));\n" ...
"  XY = [(x1 - 0.12*sin(pi*s))*Ls, d*(1-s)];\n" ...
"end\n"]);

wf(fullfile(pdir, 'GI13.m'), [ ...
"function I = GI13(M, m, n, o, type)\n" ...
"  switch type\n" ...
"    case 'Gn',  Ax=P(m); Ay=Q(n); Az=speye(o);\n" ...
"    case 'Gc',  Ax=P(m); Ay=speye(n); Az=Q(o);\n" ...
"    case 'Ge',  Ax=Q(m); Ay=P(n); Az=speye(o);\n" ...
"    case 'Gcy', Ax=speye(m); Ay=P(n); Az=Q(o);\n" ...
"    case 'Gee', Ax=Q(m); Ay=speye(n); Az=P(o);\n" ...
"    case 'Gnn', Ax=speye(m); Ay=Q(n); Az=P(o);\n" ...
"    otherwise, error('GI13: unknown type %s', type);\n" ...
"  end\n" ...
"  I = kron(Az, kron(Ay, Ax)) * M;\n" ...
"end\n" ...
"function A = Q(N)\n" ...
"  A = spdiags(0.5*ones(N,2), [0 1], N, N+1);\n" ...
"end\n" ...
"function A = P(N)\n" ...
"  A = spdiags(0.5*ones(N+1,2), [-1 0], N+1, N);\n" ...
"  A(1,1)=1.5; A(1,2)=-0.5; A(N+1,N)=1.5; A(N+1,N-1)=-0.5;\n" ...
"end\n"]);

% ---------------------------------------------------------------- setup
addpath(MOLE_SRC);
addpath(sdir);
global BT BB; BT = 2.0; BB = 5.5;

printf('\n======================================================\n');
printf('MOLE      : %s\n', MOLE_SRC);
printf('GI13 in use: %s\n', which('GI13'));
printf('Octave    : %s   (%s)\n', version(), computer());
printf('======================================================\n');

% ================================================================ CHECK 1
printf('\n===== 1. GI13 index map (no numerics) =====\n');
m=4; n=3; o=3;
[I,J,K] = ndgrid(1:m, 1:n+1, 1:o);
v = I(:) + 100*J(:) + 10000*K(:);
W = reshape(GI13(speye(m*(n+1)*o), m, n, o, 'Gn')*v, m+1, n, o);
bad = 0;
for k = 1:o
  for j = 1:n
    for i = 1:m+1
      val = W(i,j,k);
      if val == 0, bad = bad + 1; continue; end
      kk = floor(val/10000); jj = floor((val-10000*kk)/100);
      ii = val - 10000*kk - 100*jj;
      if kk ~= k
        if bad < 4
          printf('  out(%d,%d,%d) <- src(%d,%d,%d)  *** wrong zeta-plane\n', i,j,k,ii,jj,kk);
        end
        bad = bad + 1;
      end
    end
  end
end
printf('  %d of %d entries read from the wrong zeta-plane', bad, (m+1)*n*o);
printf('   (stock: 15, patched: 0)\n');

% ================================================================ CHECK 2
printf('\n===== 2. grad3DCurv convergence =====\n');
  function r = sweep()
    r = [];
    for N = [13 21 33 49]
      amp = 0.10; k = 2;
      [a,b,c] = meshgrid(linspace(0,1,N), linspace(0,1,N), linspace(0,1,N));
      X = a+amp*sin(2*pi*b); Y = b+amp*sin(2*pi*a);
      Z = c+0.5*amp*sin(2*pi*a).*sin(2*pi*b);
      s = linspace(0,1,N); s = [0, s(1:end-1)+0.5/(N-1), 1];
      [p,q,t] = meshgrid(s,s,s);
      x = p+amp*sin(2*pi*q); y = q+amp*sin(2*pi*p);
      z = t+0.5*amp*sin(2*pi*p).*sin(2*pi*q);
      f = reshape(permute(x.^2+y.^2+z.^2,[2 1 3]),[],1);
      G = grad3DCurv(k,X,Y,Z); T = G*f;
      cx = linspace(0,1,N); cx = cx(1:end-1)+0.5/(N-1);
      [aa,bb2,~] = meshgrid(linspace(0,1,N), cx, cx);
      xf = aa+amp*sin(2*pi*bb2);
      gx = T(1:N*(N-1)*(N-1));
      eg = sqrt(mean((gx - reshape(permute(2*xf,[2 1 3]),[],1)).^2));
      L = div3DCurv(k,X,Y,Z)*G; vv = reshape(L*f, N+1,N+1,N+1);
      d = vv(3:end-2,3:end-2,3:end-2)(:) - 6;
      r = [r; N, eg, sqrt(mean(d.^2))];
    end
  end
  function show(r, tag)
    printf('\n  %s\n', tag);
    printf('  %5s %13s %7s %13s %7s\n','n','grad rms','order','L=D*G rms','order');
    for i = 1:rows(r)
      if i == 1, og='   -'; ol='   -';
      else
        rr = log((r(i,1)-1)/(r(i-1,1)-1));
        og = sprintf('%5.2f', log(r(i-1,2)/r(i,2))/rr);
        ol = sprintf('%5.2f', log(r(i-1,3)/r(i,3))/rr);
      end
      printf('  %5d %13.4e %7s %13.4e %7s\n', r(i,1), r(i,2), og, r(i,3), ol);
    end
  end
show(sweep(), 'STOCK GI13');
addpath(pdir);
printf('\n  now using: %s\n', which('GI13'));
show(sweep(), 'PATCHED GI13');
printf('\n  Cartesian regression (must stay at roundoff):\n');
for N = [11 17 25]
  k = 2;
  [X,Y,Z] = meshgrid(linspace(0,1,N),linspace(0,1,N),linspace(0,1,N));
  s = linspace(0,1,N); s = [0, s(1:end-1)+0.5/(N-1), 1];
  [x,y,z] = meshgrid(s,s,s);
  f = reshape(permute(x.^2+y.^2+z.^2,[2 1 3]),[],1);
  L = div3DCurv(k,X,Y,Z)*grad3DCurv(k,X,Y,Z);
  vv = reshape(L*f, N+1,N+1,N+1);
  d = vv(3:end-2,3:end-2,3:end-2)(:) - 6;
  printf('    n=%2d  rms %.3e  max %.3e\n', N, sqrt(mean(d.^2)), max(abs(d)));
end

% ================================================================ CHECK 3
printf('\n===== 3. gridGen vs jacobian2D orientation =====\n');
here = pwd; cd(gdir);
[Xg, Zg] = gridGen('TFI', 'seamount', 99, 38, false);
cd(here);
printf('  gridGen returns  : %s   (xi-nodes as ROWS)\n', mat2str(size(Xg)));
printf('  grad2DCurv wants : rows = eta, cols = xi\n');
for c = {{'as returned', Xg, Zg}, {'transposed ', Xg.', Zg.'}}
  Jj = jacobian2D(2, c{1}{2}, c{1}{3});
  if all(Jj>0), s='ALL +'; elseif all(Jj<0), s='ALL -';
  else s=sprintf('MIXED(%d)', sum(Jj<0)); end
  printf('  %s : jacobian2D -> %-8s |J| %.4g .. %.4g\n', ...
         c{1}{1}, s, min(abs(Jj)), max(abs(Jj)));
end

% ================================================================ CHECK 4
printf('\n===== 4. TTM grid folding (thesis Eq. 2.24 seamount) =====\n');
printf('  %6s %20s %14s\n', 'iters', 'Jacobian', 'min J');
for it = [5 10 25 50 100 300]
  cd(gdir); [X,Z] = gridGen('TTM','seamount',99,38,false,it); cd(here);
  [xh,xx] = gradient(X); [zh,zx] = gradient(Z);
  Jj = xx.*zh - xh.*zx;
  if all(Jj(:)>0), s='ALL +';
  else s = sprintf('MIXED (%d neg)', sum(Jj(:)<0)); end
  printf('  %6d %20s %14.4g\n', it, s, min(Jj(:)));
end
cd(gdir); [X,Z] = gridGen('TFI','seamount',99,38,false); cd(here);
[xh,xx] = gradient(X); [zh,zx] = gradient(Z);
Jj = xx.*zh - xh.*zx;
if all(Jj(:)>0), printf('  TFI, same curves : ALL +  (usable)\n');
else printf('  TFI : MIXED (%d neg)\n', sum(Jj(:)<0)); end

printf('\n===== done =====\n\n');

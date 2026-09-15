% test_GI13_indexmap.m
%
% Demonstrates that GI13 reads from the wrong (eta, zeta) plane.
% No numerics involved -- this is purely about which source entry each
% output entry is copied from.
%
% Set MOLE_SRC below and run.  Works in Octave and MATLAB.

MOLE_SRC = '../src/matlab_octave';     % <-- adjust
addpath(MOLE_SRC);

m = 4; n = 3; o = 3;                   % CELL counts (small, so output is readable)

% Build a field on the y-faces whose value encodes its own index.
% y-faces are (m, n+1, o), flattened xi-fastest, then eta, then zeta.
[I, J, K] = ndgrid(1:m, 1:n+1, 1:o);
v = I(:) + 100*J(:) + 10000*K(:);

% Extract the raw operator by feeding GI13 an identity.
Iop = GI13(speye(m*(n+1)*o), m, n, o, 'Gn');

% Apply it.  Output lives on the x-faces, (m+1, n, o).
W = reshape(Iop*v, m+1, n, o);

fprintf('GI13(..., ''Gn'') :  y-faces (%d,%d,%d)  ->  x-faces (%d,%d,%d)\n\n', ...
        m, n+1, o, m+1, n, o);
fprintf('An output entry at (i,j,k) must be built only from source entries\n');
fprintf('on the SAME zeta-plane, i.e. k'' == k.\n\n');

bad = 0;
for k = 1:o
  for j = 1:n
    for i = 1:m+1
      val = W(i,j,k);
      if val == 0
        fprintf('  out(%d,%d,%d)  <-  ZERO\n', i, j, k);
        bad = bad + 1;
        continue
      end
      kk = floor(val/10000);
      jj = floor((val - 10000*kk)/100);
      ii = val - 10000*kk - 100*jj;
      if kk ~= k
        fprintf('  out(%d,%d,%d)  <-  src(%d,%d,%d)   *** k''=%d, expected k=%d\n', ...
                i, j, k, ii, jj, kk, kk, k);
        bad = bad + 1;
      end
    end
  end
end

fprintf('\n%d of %d output entries read from the wrong zeta-plane.\n', ...
        bad, (m+1)*n*o);
fprintf('\nCause: GI13 builds kron(speye(n*o), I1), walking the source in n*o\n');
fprintf('blocks, but the source has (n+1)*o blocks.  The trailing zero-column\n');
fprintf('pad is appended at the end rather than one plane being skipped per\n');
fprintf('zeta-level, so the map drifts by one eta-plane for every zeta-level.\n');
fprintf('The corrupted fraction grows with o.\n');

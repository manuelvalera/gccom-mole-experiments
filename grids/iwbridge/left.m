function XY = left(s)
% Lateral boundary.  The bulge that makes this a genuinely curvilinear grid
% is scaled by the DEPTH, not by Lx.  Using 0.05*Lx made the side walls
% deform without limit as the domain lengthened (150 m at Lx=3 km, 750 m at
% 15 km) while the depth stayed at 1000 m -- the resulting distortion goes
% marginal and the run becomes machine-dependent (diverged on one machine,
% collapsed to 1e-5 on another, same parameters).
    [Lx, D0, ab, Lb] = geom();
    bulge = 0.15*D0;
    XY = [-Lx/2 + bulge*sin(pi*s), -D0*(1-s)];
end

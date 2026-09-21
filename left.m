function XY = left(s)
% Lateral boundary.  The bulge that makes this a genuinely curvilinear grid
% is scaled by the DEPTH, not by Lx.  Using 0.05*Lx made the side walls
% deform without limit as the domain lengthened (150 m at Lx=3 km, 750 m at
% 15 km) while the depth stayed at 1000 m -- the resulting distortion goes
% marginal and the run becomes machine-dependent (diverged on one machine,
% collapsed to 1e-5 on another, same parameters).
%
% BULGE is a global so the driver can set it without rewriting this file.
% Hand-editing it is how the file got left at bulge=0 once, which silently
% turned several runs into a different experiment.  Default 0.15 reproduces
% the original behaviour exactly.
    global BULGE
    if isempty(BULGE), BULGE = 0.15; end
    [Lx, D0, ab, Lb] = geom();
    bulge = BULGE*D0;
    XY = [ -Lx/2 + bulge*sin(pi*s), -D0*(1-s)];
end

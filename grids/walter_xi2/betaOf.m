function bb = betaOf(s, which)
% thesis 5.1: beta varies with position, and differs top vs bottom, so the
% xi-lines are NOT vertical -- a fully curvilinear grid, not sigma.
    global BT BB
    if isempty(BT), BT=1.6; end
    if isempty(BB), BB=3.2; end
    if strcmp(which,'bot'), bb = BB*(1 - 0.30*sin(pi*s));
    else,                   bb = BT*(1 + 0.45*sin(pi*s));
    end
end

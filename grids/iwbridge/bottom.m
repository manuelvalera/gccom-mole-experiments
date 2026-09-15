function XY = bottom(s)
% Garcia et al. 2019 sec 3.3 ridge: D = D0 - a_b exp(-x^2/(2 L_b^2))
    [Lx, D0, ab, Lb] = geom();
    x = -Lx/2 + Lx*stretch(s,1,betaOf(s,'bot'));
    XY = [x, -(D0 - ab*exp(-x^2/(2*Lb^2)))];
end

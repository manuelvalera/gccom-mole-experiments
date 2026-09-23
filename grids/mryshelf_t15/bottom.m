function XY = bottom(s)
    [Lx, D0, ab, Lb] = geom();
    x = -Lx/2 + Lx*stretch(s, 1, betaOf(s, 'bot'));
    XY = [x, -mry_depth(x)];
end

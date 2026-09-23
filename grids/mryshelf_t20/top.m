function XY = top(s)
    [Lx, D0, ab, Lb] = geom();
    x = -Lx/2 + Lx*stretch(s, 1, betaOf(s, 'top'));
    XY = [x, 0];
end

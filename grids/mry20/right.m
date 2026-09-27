function XY = right(s)
    [Lx, D0, ab, Lb] = geom();
    global BULGE
    if isempty(BULGE)
        BULGE = 0.0;
    end
    d = mry_depth(Lx/2);
    XY = [Lx/2 - BULGE*D0*sin(pi*s), -d*(1 - s)];
end

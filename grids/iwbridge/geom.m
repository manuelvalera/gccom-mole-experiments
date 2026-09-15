function [Lx, D0, ab, Lb] = geom()
% Domain geometry, overridable by the driver.  gridGen resolves boundary
% curves by name from the working directory, so the driver writes geom_over.m
% next to this file when --Lx / --D0 are given.
    if exist('geom_over') == 2
        [Lx, D0, ab, Lb] = geom_over();
    else
        Lx = 3000; D0 = 1000; ab = 20; Lb = 30;
    end
end

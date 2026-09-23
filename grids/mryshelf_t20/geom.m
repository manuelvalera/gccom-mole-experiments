function [Lx, D0, ab, Lb] = geom()
% Monterey shelf transect. ab and Lb are unused here (the bed comes
% from measured bathymetry) but are kept so the signature matches the
% other grids and the driver's geom_over mechanism still works.
    if exist('geom_over') == 2
        [Lx, D0, ab, Lb] = geom_over();
    else
        Lx = 9409.3;
        D0 = 87.78;
        ab = 0;
        Lb = 1;
    end
end

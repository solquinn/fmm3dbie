S = geometries.disk([],[],[3 3 3],8);

thetas = (0:23)/12*pi;
targinfo = []; targinfo.r = [cos(thetas); sin(thetas)];
targinfo.n = targinfo.r;

A = lap2d.v2b_neu(S,targinfo,1e-9);
A = reshape(A,[S.npts length(targinfo.r)]).';

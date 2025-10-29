S = geometries.disk([],[],[5 5 5],8);

targinfo = []; targinfo.r = [-1;1;0];
targinfo.n = targinfo.r;

dens = test_fn(S.r(1,:),S.r(2,:)).';

A = lap2d.evalmat(S,targinfo,1e-12);
A = reshape(A,[S.npts size(targinfo.r,2)]).';

val = integral2(@(r,t) -1/(2*pi)*r.*log(sqrt((r.*cos(t) - targinfo.r(1)).^2 + (r.*sin(t) - targinfo.r(2)).^2)).*test_fn(r.*cos(t),r.*sin(t)), 0,1,0,2*pi,"AbsTol",1e-12,"RelTol",1e-12);

val2 = A*dens;

abs(val - val2) / abs(val)

function val = test_fn(x,y)

    val = exp( - 100*x.^2 - 50*y.^4 );

end
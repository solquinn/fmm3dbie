S = geometries.disk([],[],[3 3 3],8);

if false

thetas = pi/24;

% thetas = 0;
targinfo = []; targinfo.r = [cos(thetas); sin(thetas); 0];
targinfo.n = targinfo.r;

tic;
Av2b = lap2d.v2b_neu(S,targinfo,1e-8);
toc; 

h = 1e-6;
targinfo2 = []; targinfo2.r = targinfo.r + h*(0:2).*targinfo.n;

tic;
Aeval = lap2d.evalmat(S,targinfo2,1e-8);
toc;

dens = test_fn(S.r(1,:),S.r(2,:)).';

val = Av2b*dens;

d1 = [-3/2	2	-1/2]/h;
val2 = d1*(Aeval*dens);

abs(val - val2) / abs(val)

end

function val = test_fn(x,y)

    val = exp( - 100*x.^2 - 50*y.^4 );

end
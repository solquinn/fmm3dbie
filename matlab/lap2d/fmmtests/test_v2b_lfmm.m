S = geometries.disk([],[],[4 4 4],8);

thetas = (0:12)*pi/24;

targinfo = []; targinfo.r = [cos(thetas); sin(thetas); 0*thetas];
targinfo.n = targinfo.r;

tic;
Av2b = lap2d.v2b_neu(S,targinfo,1e-12);
toc; 

rs = S.r;
test_fn = exp(-5*((rs(1,:) - 0.3).^2+(rs(2,:)+1/pi).^2));
test_fn = test_fn.';

sol1 = Av2b*test_fn;

% lap2d_s = kernel('l','sp');
% A_native = lap2d_s.eval(S_over,targinfo).*S_over.wts(:).'; 
% sol2 = A_native*test_fn_over + Av2b_cor*test_fn;

tic;
[Av2b_cor,nover] = lap2d.get_quad_cor_v2b_neu(S, targinfo, 1e-12);
toc;

sol2 = lap2d.apply_v2b_neu(S,targinfo,test_fn,Av2b_cor,nover,1e-12);

vecnorm(sol1 - sol2) / vecnorm(sol1)

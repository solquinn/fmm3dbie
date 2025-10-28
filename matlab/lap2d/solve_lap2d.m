S = geometries.disk([],[],[3 3 3],8);

figure(1); clf
plot(S,rand(S.npatches,1))

%%
A = lap2d.slp_matgen(S,1e-9);
A = reshape(A,[S.npts,S.npts]).';

%%

V = eval_gauss(S.r);
rhs = eval_gauss(S.r);

lhs = -eye(S.npts) + V.*A;

sigma = lhs \ rhs;

t = tiledlayout(1,4);
nexttile
plot(S,rhs)
colorbar
view(0,90)
title('f')

nexttile
plot(S,V)
colorbar
view(0,90)
title('V')

nexttile
plot(S,sigma)
colorbar
view(0,90)
title('\sigma')

u = A* sigma;
resid = get_surface_laplacian(S,u) + V.*u - rhs;

nexttile 
scatter(S.r(1,:),S.r(2,:),8,log10(abs(resid)))
colorbar
title('log_{10} residual')

function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end
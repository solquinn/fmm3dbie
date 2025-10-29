S = geometries.disk([],[],[3 3 3],8);

figure(1); clf
plot(S,rand(S.npatches,1))

A = lap2d.slp_matgen(S,1e-9);
A = reshape(A,[S.npts,S.npts]).';

tic;
[v2v_cor,nover] = lap2d.get_quad_cor_sub(S, 1e-12);
toc;

v2v_apply = @(mu) apply_lap2d_v2v(S,mu,v2v_cor,nover);

V = eval_gauss(S.r);
rhs = eval_gauss(S.r);

lhsmat = -eye(S.npts) + V.*A;
lhs = @(mu) -mu + V.*v2v_apply(mu);

mu = gmres(@(mu) lhs(mu),rhs,[],1e-10,2000);
u = v2v_apply(mu);

resid = abs(get_surface_laplacian(S,u) + V.*u - rhs) / max(abs(mu));


t = tiledlayout(1,2);
nexttile
plot(S,u)
colorbar
view(0,90)
title('u')

nexttile 
scatter(S.r(1,:),S.r(2,:),8,log10(abs(resid)))
axis square
colorbar
title('log_{10} residual')

function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end
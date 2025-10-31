S = geometries.disk([],[],[4 4 4],12);

figure(1); clf
plot(S,rand(S.npatches,1))

tic;
[v2v_cor,nover] = lap2d.get_quad_cor_sub(S, 1e-12);
toc;

v2v_apply = @(mu) lap2d.apply_v2v(S,mu,v2v_cor,nover,1e-12);

V = eval_V(S.r);
rhs = eval_rhs(S.r);

lhs = @(mu) -mu + V.*v2v_apply(mu);

mu = gmres(@(mu) lhs(mu),V.*rhs,[],1e-10,2000);
u = v2v_apply(mu);

resid = abs(get_surface_laplacian(S,u) + V.*u - V.*rhs);

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


function V = eval_V(targ)

V = exp( - 10*(targ(1,:)).^2 - 10*(targ(2,:)).^2 );
V = V(:);

end

function rhs = eval_rhs(targ)

rhs = exp( - 10*(targ(1,:)-0.25).^2 - 10*(targ(2,:)-0.25).^2 );
rhs = rhs(:);

end
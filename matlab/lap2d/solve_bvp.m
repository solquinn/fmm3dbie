% genpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/')
% addpath('/Users/squinn/chunkie/chunkie/')

S = geometries.disk([],[],[4 4 4],8);

% chnkr = chunkerfunc(@(t) starfish(t,5,0,[0,0],0,1));
nch = 4*4;
cparams.ta = pi/nch;
cparams.tb = 2*pi + cparams.ta;
chnkr = chunkerfuncuni(@(t) ellipse(t),nch,cparams);
chnkr = sort(chnkr);

figure(1); clf
plot(S,rand(S.npatches,1))
hold on
plot(chnkr,'x-')
view(0,90)

%%
V = eval_gauss(S.r);
%rhs = eval_gauss(xs);

%%
% A = lap2d.slp_matgen(S,1e-9);
% lhs_11 = -eye(S.npts) + V.*A;

tic;
[Av2v_cor,nover] = lap2d.get_quad_cor_sub(S, 1e-12);
toc;

v2v_apply = @(mu) lap2d.apply_v2v(S,mu,Av2v_cor,nover,1e-12);
lhs_11 = @(mu) - mu + V.*v2v_apply(mu);

%%

%fkern = chnk.lap2d.kern;
fkern = @(s,t) chnk.lap2d.kern(s,t,'s');
lhs_12 = V.*chunkerkernevalmat(chnkr,fkern,S.r(1:2,:));

%%
targinfo=[];
targinfo.r = [chnkr.r(:,:);0*chnkr.r(1,:)];
targinfo.n = [chnkr.n(:,:);0*chnkr.n(1,:)];
% lhs_21 = lap2d.v2b_neu(S,targinfo,1e-8);

tic;
[Av2b_cor,nover] = lap2d.get_quad_cor_v2b_neu(S, targinfo, 1e-12);
toc;

v2b_apply = @(mu) lap2d.apply_v2b_neu(S,targinfo,mu,Av2b_cor,nover,1e-12);
lhs_21 = @(mu) v2b_apply(mu);

%%
fkern_prime = @(s,t) chnk.lap2d.kern(s,t,'sprime');
lhs_22 = 0.5*eye(chnkr.npt)+chunkermat(chnkr,fkern_prime); %0.5*eye(n)... -

%%

% lhs = [lhs_11, lhs_12; lhs_21, lhs_22];
lhs = @(x) [lhs_11(x(1:S.npts)) + lhs_12*x(S.npts+1:end); lhs_21(x(1:S.npts)) + lhs_22*x(S.npts+1:end)];

%rhs_1 = eval_gauss(xs);
%rhs_2 = zeros(n,1);

% analytic solution : u = sin(x)sin(y)
rhs_1 = (-2+V(:).').*sin(S.r(1,:)).*sin(S.r(2,:)) ; 
rhs_2 = cos(chnkr.r(1,:)).*sin(chnkr.r(2,:)).*chnkr.n(1,:) ...
          + sin(chnkr.r(1,:)).*cos(chnkr.r(2,:)).*chnkr.n(2,:) ; % sin(xs(1,:))

rhs = [rhs_1 rhs_2].';

% dens = lhs\rhs;
dens = gmres(lhs,rhs,[],1e-10,2000);
mu = dens(1:S.npts);
rho = dens(S.npts+1:end);

fkern_s = kernel('l','s');
u = v2v_apply(mu) + chunkerkerneval(chnkr,fkern_s,rho,S.r(1:2,:));

ref_u = sin(S.r(1,:)).*sin(S.r(2,:));
err = abs(u - ref_u(:)) / max(abs(u));

figure(1); clf
scatter(S.r(1,:),S.r(2,:),8,log10(err)); colorbar

fprintf('relative L_2/L_2 error: %5.5e\n', vecnorm(err .* S.wts(:)) / vecnorm(u .* S.wts(:)) );

return

%%
% targvol = [];
% targvol.r = xs;
% B = chunkerkernevalmat(chnkr,fkern,targvol);
% B = chunkerkernevalmat(chnkr,fkern,xs);
% u_sol = A*dens_1 + B*dens_2;

%%
if false
%%
%test PDE
resid = get_surface_laplacian(S,u_sol) + V.*u_sol - rhs_1;
scatter(S.r(1,:),S.r(2,:),8,log10(abs(resid)))
axis square
colorbar
title('log_{10} residual')

%%
%test BCs
h = 1e-3;
bdry_idx = [56]; %just one index for now
bdry_pts = rs(:,bdry_idx);
bdry_normals = targinfo.n(1:2,bdry_idx);
targinfo2 = []; targinfo2.r = [1,1,1].*bdry_pts + h*(0:2).*bdry_normals;
%%
tic;
Aeval = lap2d.evalmat(S,targinfo2,1e-8);
toc;
%%
d1 = [-3/2	2	-1/2]/h;
val2 = d1*(Aeval*mu)

end

%%
%analytic solution test
u_true = 0.5*(xs(1,:).^2) - xs(2,:);
u_true = u_true.';
err = u_sol - u_true;
scatter(S.r(1,:),S.r(2,:),8,log10(abs(err)))
axis square
colorbar
title('log_{10} error')

%%
figure;
scatter(S.r(1,:),S.r(2,:),8,real(u_true))
axis square
colorbar
title('true solution')

%%
figure;
scatter(S.r(1,:),S.r(2,:),8,real(u_sol))
axis square
colorbar
title('computed solution')

%%
figure;
scatter(S.r(1,:),S.r(2,:),8,real(mu))
axis square
colorbar
title('volumetric density')

%%
figure;
scatter(S.r(1,:),S.r(2,:),8,real(rhs_1))
axis square
colorbar
title('$f$','Interpreter','latex')

%%
if false

%%
[U, singular_vals, V] = svd(lhs);

%%
singular_vals = diag(singular_vals);

%%
figure;
plot(S.r(1,:), S.r(2,:),'.');
hold on
plot(rs(1,:), rs(2,:),'x');

%%
figure;
plot(chnkr.r(1,:),real(rho))
%%
figure;
plot(S, real(mu)); colorbar
figure; plot(S, rhs_1); colorbar

%%
figure;
plot(S, real(u_sol)); colorbar

end

%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end
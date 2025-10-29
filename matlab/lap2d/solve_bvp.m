genpath('/Users/squinn/chunkie/')
%addpath('/Users/squinn/chunkie/')
addpath('/Users/squinn/chunkie/chunkie/')

S = geometries.disk([],[],[3 3 3],8);

%figure(1); clf
%plot(S,rand(S.npatches,1))

xs= S.r(1:2,:);

%%
V = eval_gauss(xs);
%rhs = eval_gauss(xs);

%%
A = lap2d.slp_matgen(S,1e-9);
A = reshape(A,[S.npts,S.npts]).'; %deleted minus sign.

%%
lhs_11 = -eye(S.npts) + V.*A;

%%
chnkr = chunkerfunc(@(t) starfish(t,5,0,[0,0],0,1));
temp=chnkr.r;
rs = reshape(temp,2,[]);

n = chnkr.k*chnkr.nch;

normals = chnkr.d;
normals = reshape(normals,2,[]);

%fkern = chnk.lap2d.kern;
fkern = @(s,t) chnk.lap2d.kern(s,t,'s');
lhs_12 = V.*chunkerkernevalmat(chnkr,fkern,xs);

%%
targinfo=[];
targinfo.r = [rs;zeros(1,n)];
targinfo.n = [normals;zeros(1,n)];
lhs_21 = lap2d.v2b_neu(S,targinfo,1e-8);

%%
fkern_prime = @(s,t) chnk.lap2d.kern(s,t,'sprime');
lhs_22 = 0.5*eye(n)+chunkermat(chnkr,fkern_prime); %0.5*eye(n)... -

%%
lhs = [lhs_11, lhs_12; lhs_21, lhs_22];

%rhs_1 = eval_gauss(xs);
%rhs_2 = zeros(n,1);
xs_norms = sqrt(xs(1,:).^2+xs(2,:).^2);
rhs_1 = 1+exp(-10*xs_norms.^2).*(0.5*xs(1,:).^2 - xs(2,:));
rhs_2 = rs(1,:).^2 - rs(2,:);
rhs_1 = rhs_1.'; rhs_2 = rhs_2.';

rhs = [rhs_1; rhs_2];

dens = lhs\rhs;
dens_1 = dens(1:S.npts);
dens_2 = dens(S.npts+1 : S.npts+n);

%%
% targvol = [];
% targvol.r = xs;
% B = chunkerkernevalmat(chnkr,fkern,targvol);
B = chunkerkernevalmat(chnkr,fkern,xs);
u_sol = A*dens_1 + B*dens_2;

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
val2 = d1*(Aeval*dens_1)

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
scatter(S.r(1,:),S.r(2,:),8,real(dens_1))
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
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end

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
plot(chnkr.r(1,:),real(dens_2))
%%
figure;
plot(S, real(dens_1)); colorbar
figure; plot(S, rhs_1); colorbar

%%
figure;
plot(S, real(u_sol)); colorbar

end
S = geometries.disk([],[],[4 4 4],8);

%figure(1); clf
%plot(S,rand(S.npatches,1))

%%

rs = S.r;
test_fn = exp(-5*(rs(1,:).^2+rs(2,:).^2));
test_fn = test_fn.';

tic;
A = lap2d.slp_matgen(S,1e-10);
A = reshape(A,[S.npts,S.npts]).';
toc;

sol1 = A*test_fn;

tic;
[v2v_cor,nover] = lap2d.get_quad_cor_sub(S, 1e-12);
toc;

[S_over,xinterp] = oversample(S,nover);
test_fn_over = xinterp*test_fn;

% lap2d_s = kernel('l','s');
% A_native = lap2d_s.eval(S_over,S).*S_over.wts(:).';
% 
% sol2 = A_native*test_fn_over + v2v_cor*test_fn;

srcinfo = []; 
srcinfo.sources = S_over.r(1:2,:);
srcinfo.charges = test_fn_over(:).'.*S_over.wts(:).';
targ = S.r(1:2,:);
U = lfmm2d(eps,srcinfo,0,targ,1);
sol2 = -1/(2*pi)*U.pottarg.';
sol2 = sol2 + v2v_cor*test_fn(:);

vecnorm(sol1 - sol2) / vecnorm(sol1)

return

%%


%%
dens = A*test_fn;

lapdens = get_surface_laplacian(S,-dens);
err = abs(test_fn - lapdens);

fprintf('relative l_inf/l_1 error interior: %5.2e\n', vecnorm(err,'Inf') / sum(abs(test_fn(:) .* S.wts(:))));

%%
figure; 
scatter(S.r(1,:),S.r(2,:),8,log10(err)); colorbar
axis square


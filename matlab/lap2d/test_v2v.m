S = geometries.disk([],[],[4 4 4],8);

%figure(1); clf
%plot(S,rand(S.npatches,1))

%%
A = lap2d.slp_matgen(S,1e-10);
A = reshape(A,[S.npts,S.npts]).';

%%
rs = S.r;
test_fn = exp(-5*(rs(1,:).^2+rs(2,:).^2));
test_fn = test_fn.';

%%
dens = A*test_fn;

lapdens = get_surface_laplacian(S,-dens);
err = abs(test_fn - lapdens);

fprintf('relative l_inf/l_1 error interior: %5.2e\n', vecnorm(err,'Inf') / sum(abs(test_fn(:) .* S.wts(:))));

%%
figure; 
scatter(S.r(1,:),S.r(2,:),8,log10(err)); colorbar
axis square


S = geometries.disk([],[],[3 3 3],8);

%figure(1); clf
%plot(S,rand(S.npatches,1))

%%
A = lap2d.slp_matgen(S,1e-9);
A = reshape(A,[S.npts,S.npts]).';

%%
rs = S.r;
test_fn = exp(-5*(rs(1,:).^2+rs(2,:).^2));
test_fn = test_fn.';

%%
dens = A*test_fn;

%%
figure; plot(S,test_fn)
figure; plot(S,dens)


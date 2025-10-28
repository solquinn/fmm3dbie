genpath('/Users/squinn/chunkie/')
%addpath('/Users/squinn/chunkie/')
addpath('/Users/squinn/chunkie/chunkie/')

S = geometries.disk([],[],[3 3 3],8);

%figure(1); clf
%plot(S,rand(S.npatches,1))

xs= S.r([1:2],:);

%%
V = eval_gauss(xs);
rhs = eval_gauss(xs);

%%
A = lap2d.slp_matgen(S,1e-9);
A = -reshape(A,[S.npts,S.npts]).'; %added minus sign.

%%
lhs_11 = eye(S.npts) + V.*A;

%%
chnkr = chunkerfunc(@(t) starfish(t,5,0.3,[0,0],0,2));
temp=chnkr.r;
rs = reshape(temp,2,[]);

%fkern = chnk.lap2d.kern;
fkern = @(s,t) chnk.lap2d.kern(s,t,'s');
lhs_12 = V.*chunkerkernevalmat(chnkr,fkern,xs);

%%
lhs_21 = 

%%
n = chnkr.k*chnkr.nch;
lhs_22 = -0.5*eye(n)+chunkermat(chnkr,fkern);

%%
function val = eval_gauss(targ)

val = exp( - 10*targ(1,:).^2 - 10*targ(2,:).^2 );
val = val(:);

end
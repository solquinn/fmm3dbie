function [val] = green(src,targ)
%HELM3D.GREEN evaluate the laplace green's function
% for the given sources and targets

[~,ns] = size(src);
[~,nt] = size(targ);

xs = repmat(src(1,:),nt,1);
ys = repmat(src(2,:),nt,1);

xt = repmat(targ(1,:).',1,ns);
yt = repmat(targ(2,:).',1,ns);

rx = xt-xs;
ry = yt-ys;

rx2 = rx.*rx;
ry2 = ry.*ry;

r2 = rx2+ry2;

r = sqrt(r2);

val = -(1/(2*pi))*log(r);

end
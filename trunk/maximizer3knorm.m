function [u,v,w,ev] = maximizer3knorm(varargin)
% MAXIMIZER3KNORM - find best rank-1 approximation to a tensor in 3-norm
% maximizer3 - maximizes sum_{i,j,k} u(i)v(j)w(k) T(i,j,k) subject to the
% constraints norm(u,p)=norm(v,p)=norm(w,p)=1;
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

T = varargin{1};
status = varargin{2};
p = varargin{3};

% check if there is  motif generalization
if isempty(T.I) 
u = sparse(T.dim(1),1);
v = sparse(T.dim(2),1);
w = sparse(T.dim(3),1);
ev = 0;
return;
end

% check if initial u, v are given, else use default
if nargin >= 5
u = varargin{4};
v = varargin{5};
else
% initialize u, v to pair with highest number of common motifs
[I,J,Av]=find(tensmat(T,ones(T.dim(3),1),3));
[amax,imax] = max(Av);
u = sparse(T.dim(1),1);
v = sparse(T.dim(1),1);
u(I(imax)) = 1;
v(J(imax)) = 1;
end

% check if initial w is given, else use default
if nargin == 6
w = varargin{6};
else
% initialize w by Euler-Lagrange or symmetry
w = tensmat(T,u,v,[1 2]).^(1/(p-1));
w = w/norm(w,p);
end

% symmetrize if necessary
switch status
case 1
u = (u+v+w)/norm(u+v+w,p);
v = u;
w = u;
case 2
v = (v+w)/norm(v+w,p);
w = v;
end

% initialize convergence tester
diff = 1;
step = 0;
% loop until convergence or maximum number of steps
tol = 1e-5;
maxstep = 500;
while diff > tol && step<maxstep
step = step+1;
% update vectors
unew = tensmat(T,v,w,[2 3]).^(1/(p-1));
unew = unew/norm(unew,p);
vnew = tensmat(T,unew,w,[1 3]).^(1/(p-1));
vnew = vnew/norm(vnew,p);
wnew = tensmat(T,unew,vnew,[1 2]).^(1/(p-1));
wnew = wnew/norm(wnew,p);
% set eigenvalue
ev = tensmat(T,unew,vnew,wnew,[1 2 3]);
% update convergence tester
diff = max(max(norm(unew-u),norm(vnew-v)),norm(wnew-w));
% reset
u = unew;
v = vnew;
w = wnew;
% symmetrize if necessary
switch status
case 1
u = (u+v+w)/norm(u+v+w,p);
v = u;
w = u;
case 2
v = (v+w)/norm(v+w,p);
w = v;
end
end
if step==maxstep
warning('tomic:maximizer3', 'Maximizer3: maximum number of steps reached, diff = %1.3f', diff);
end

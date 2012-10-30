function [y,o] = unifapproxknorm(varargin)
% UNIFAPPROX - find best uniform approximation to x in k-norm sense
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

x = varargin{1};
k = varargin{2};
x = x./norm(x,k); % normalize if necessary
tol = 1e-5; % set to the same value as used in maximizer functions

% check if x is nonnegative
if sum(x<0)>0
error('Input vector should be nonnegative')
end

% sort entries of x in descending order
[xt,t] = sort(full(x),'descend');
xt(xt<=tol) = 0;
xt = sparse(xt);

% create all normalized index vectors, column l corresponds to taking first
% l largest entries of x
M = triu(repmat(spones(xt),[1,sum(xt>tol)]))./...
repmat((1:length(find(xt>tol))).^(1/k),[length(xt),1]);

% compute k-norm distance with normalized index vectors
sc = zeros(1,size(M,2));
for l=1:length(sc)
sc(l) = norm(xt-M(:,l),k); 
end
plot(sc(1:10))
% locate minimum distance
[o,pos] = min(sc);
% rescue mode for rare cases: take n-th best choice
if nargin==3
[sc2,tt] = sort(full(sc));
o = sc2(varargin{3});
pos = tt(varargin{3});
end

% get corresponding vector and return to original index ordering
y = M(:,pos);
[t1,t2] = sort(t); % t2 = inverse permutation of t
y = spones(y(t2));


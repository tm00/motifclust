function t = tenssym(t1)
% TENSSYM - symmetrize a tensor
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

ndim = size(t1.dim,2);
if (ndim>=8)
error('Symmetrization is only practical when the number of dimensions is less than 8')
end

% generate all permutations
p = perms(1:ndim);

% generate all index combinations and corresponding values
I = reshape(t1.I(:,p), [size(t1.I,1)*size(p,1), size(t1.I,2)]);
V = repmat(t1.V,[size(p,1), 1]);

% there might be duplicate entries
[I, m] = unique(I,'rows');
V = V(m);

% create output tensor
t.I = I;
t.V = V;
t.dim = t1.dim;

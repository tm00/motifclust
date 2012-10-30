function t = tensmin(t1,t2)
% TENSMIN - difference of two tensors
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

if t1.dim ~= t2.dim
error('Two tensors need to have same dimensions');
end

[I, m, n] = unique([t1.I; t2.I],'rows');
V = accumarray(n, [t1.V; -t2.V]);
I = I(V>0,:);
V = V(V>0);
t.I = I;
t.V = V;
t.dim = t1.dim;

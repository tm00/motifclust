function S = tenstrunc(s1,s2,s3,T)
% TENSTRUNC - truncate tensor
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

% create index array for kron(kron(s3,s2),s1)
u1 = find(s1);
u2 = find(s2);
u3 = find(s3);

n1 = length(u1);
n2 = length(u2);
n3 = length(u3);

I1 = repmat(u1', [n2*n3,1]);
I1 = reshape(I1, [n1*n2*n3,1]);

I2 = repmat(u2', [n3,1]);
I2 = reshape(I2, [n2*n3,1]);
I2 = repmat(I2, [n1,1]);

I3 = repmat(u3, [n1*n2,1]);

% find overlap and truncate
[I, i1, i2] = intersect2([I1, I2, I3], T.I, 'rows');
V = T.V(i2);
S.I = I;
S.V = V;
S.dim = T.dim;

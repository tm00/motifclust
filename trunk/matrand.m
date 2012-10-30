function s = matrand(a)
% MATRAND - randomize matrix
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
%  WRITTEN BY
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de

% check diagonal
if ~isempty(find(diag(a))) & size(a,1)==size(a,2)
    warning('Matrand assumes empty diagonal, removing entries on diagonal.')
    a = a - diag(diag(a));
end

% test symmetry
issym = size(a,1)==size(a,2) && isempty(find(a-a'));

% convert matrix to 2D tensor
[I,J,V] = find(a);
T.I = [I,J];
T.V = V;
T.dim = size(a);

% randomize tensor
if issym
    S = tensrand(T,[1 2]);
else
    S = tensrand(T);
end

% convert back to matrix
s = sparse(S.I(:,1),S.I(:,2),S.V,size(a,1),size(a,2));

if issym
    s = 0.5*(s+s');
end

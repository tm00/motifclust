function A = edgeswap(A,i1,j1,i2,j2)
% EDGESWAP - swap two edges in a randomization procedure
%
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

issym = size(A,1)==size(A,2) && isempty(find(A-A'));
if length(unique([i1,j1,i2,j2]))==4 && A(i1,j1)>0 && A(i2,j2)>0 && A(i1,j2)==0 && A(i2,j1)==0
    A(i1,j1) = 0;
    A(i1,j2) = 1;
    A(i2,j1) = 1;
    A(i2,j2) = 0;
    if issym
        A(j1,i1) = 0;
        A(j2,i1) = 1;
        A(j1,i2) = 1;
        A(j2,i2) = 0;
    end
end

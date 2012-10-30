function oc = aggrscore(T)
% AGGRSCORE - aggregation score of network motif tensor
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de

% number of nodes in each dimension
nd = zeros(size(T.I,2),1);
for k=1:length(nd)
    nd(k) = length(unique(T.I(:,k)));
end

% overlap coefficient
oc = size(T.I,1)/prod(nd.^0.5);

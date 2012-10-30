function [pval, Z, sreal, srand] = motif3aggregation(T, Trand)
% MOTIF3AGGREGATION - computes the significance for the aggregation of a network motif
% MOTIF3AGGREGATION computes the p-value and Z-score for the aggregation of
% a certain network motif using the overlap coefficient of the motif
%
% USAGE:
%   [pval, Z, Hreal, Hrand] = motif3aggregation(T, Trand)
%
% PARAMETERS:
%   - T : motif tensor for a given input motif
%   - Trand : cell array of motif tensors obtained by randomizing the
%     appropriate networks while keeping the motif count constant
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
%  WRITTEN BY
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de

nR = length(Trand);
srand = zeros(nR,1);
for k=1:nR
    srand(k) = aggrscore(Trand{k});
end
sreal = aggrscore(T);
Z = (sreal-mean(srand))./std(srand);
pval = sum(srand>=sreal)/length(srand);
%Zrand = (srand - repmat(mean(srand),[nR,1]))./repmat(std(srand),[nR,1]);

function [nm, sc] = motif3cluststat(varargin)
% MOTIF3CLUSTSTAT - get some statistics from motif3clust output
% MOTIF3CLUSTSTAT computes the number of motifs and the score of each
% cluster from the output of motif3clust.
%
% USAGE:
%   - [nm,sc] = motif3cluststat(x,(y,z),S)
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

switch nargin
case 2
x = varargin{1};
S = varargin{2};
nm = zeros(length(S),1);
sc = zeros(length(S),1);
for l=1:length(S)
nm(l) = length(S{l}.V)/6;
nx = sum(x(:,l)>0);
sc(l) = nm(l)/(nx^1.5);
end
case 3
x = varargin{1};
y = varargin{2};
S = varargin{3};
nm = zeros(length(S),1);
sc = zeros(length(S),1);
for l=1:length(S)
nm(l) = length(S{l}.V)/2;
nx = sum(x(:,l)>0);
ny = sum(y(:,l)>0);
sc(l) = nm(l)/(sqrt(nx)*ny);
end
case 4
x = varargin{1};
y = varargin{2};
z = varargin{3};
S = varargin{4};
nm = zeros(length(S),1);
sc = zeros(length(S),1);
for l=1:length(S)
nm(l) = length(S{l}.V);
nx = sum(x(:,l)>0);
ny = sum(y(:,l)>0);
nz = sum(z(:,l)>0);
sc(l) = nm(l)/(sqrt(nx*ny*nz));
end
end

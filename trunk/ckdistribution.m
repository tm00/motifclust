function [c,p,n,d] = ckdistribution(a)
% CKDISTRIBUTION - average degree and cluster coefficient distribution
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de

% we need symmetric edge matrix (0/1's only)
issym = size(a,1)==size(a,2) && isempty(find(a-a'));
if ~issym
    a = spones(a+a');
else
    a = spones(a);
end

[a,cc,dd] = clustcoeff(a);

cc = full(cc(dd>0));
dd = full(dd(dd>0));

% unique degree values
c = accumarray(dd,cc);
n = accumarray(dd,ones(size(dd)));

d = find(n>0);
n = n(d);
c = c(d)./n;
p = n./sum(n);

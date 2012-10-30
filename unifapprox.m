function [y,o] = unifapprox(varargin)
% UNIFAPPROX - find best uniform approximation to x in 2-norm sense
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
x = x./norm(x); % normalize if necessary
tol = 1e-5; % set to the same value as used in maximizer functions

% check if x is nonnegative
if sum(x<0)>0
error('Input vector should be nonnegative')
end

% sort entries of x in descending order
[xt,t] = sort(full(x),'descend');
xt(xt<=tol) = 0;
xt = sparse(xt);

% create all truncated vectors, column k contains vector of first k largest
% entries of x
M = triu(repmat(xt,[1,length(xt)]));

% compute l2-overlap with normalized index vectors
sc = sum(M)./sqrt(1:length(xt));

% locate maximum overlap
[o,pos] = max(sc);
% rescue mode for rare cases: take k-th best choice
if nargin==2
[sc2,tt] = sort(full(sc),'descend');
o = sc2(varargin{2});
pos = tt(varargin{2});
end

% get corresponding vector and return to original index ordering
y = M(:,pos);
[t1,t2] = sort(t); % t2 = inverse permitation of t
y = y(t2);
y = y./norm(y);


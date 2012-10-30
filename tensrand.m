function S = tensrand(varargin)
% TENSRAND - randomize 3D tensor
% Tensrand generates a random tensor with the same one-mode degree
% distributions as an input tensor T.
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

% read input
T = varargin{1};
sym = cell(nargin-1,1);
for k=1:nargin-1
sym{k} = varargin{k+1}; % each argument is a pair of dimensions with symmetry
end

% for symmetric dimensions, the input contains each pair twice; for the 
% randomization we keep each only once
for k=1:length(sym)
I = T.I(:,[sym{k}(1),sym{k}(2)]);
T.I = T.I(I(:,2)>I(:,1),:);
T.V = T.V(I(:,2)>I(:,1));
end

% number of motifs and dimensions
N = length(T.V);
ndim = length(T.dim);

% randomize
S.I = zeros(size(T.I));
for k=1:ndim
S.I(:,k) = T.I(randperm(N),k);
end
% copy dimensions etc.
S.V = T.V;
S.dim = T.dim;

% each motif should consist of distinct nodes
[pdI,pdJ] = find(triu(ones(ndim),1));
for r=1:length(pdI) % check all dimension pairs
S = correctcolumns(S,[pdI(r) pdJ(r)]);
end

nstep = 100;
for n=1:nstep
if length(unique(S.I,'rows'))==N
break;
end
% S should not contain duplicate entries
S = correctrows(S);

% for symmetric dimensions, we need to correct
for k=1:length(sym)
S = correctsym(S, sym{k});
end
end
if n==nstep
warning('Tensrand: maximum number of iterations reached.');
end

% S should not contain duplicate entries
%S = correctrows(S);

% re-symmetrize symmetric dimensions
for k=1:length(sym)
I2 = S.I;
I2(:,[sym{k}(1),sym{k}(2)]) = S.I(:,[sym{k}(2),sym{k}(1)]);
S.I = [S.I; I2];
S.V = [S.V; S.V];
end

function S = correctsym(S, c)
% correct symmetry such that entries in c(1) are always smaller than in
% c(2)
I = S.I(:,c);
ix = I(:,1)>I(:,2);
% switch columns where necessary 
S.I(ix,c(2)) = I(ix,1);
S.I(ix,c(1)) = I(ix,2);

function S = correctcolumns(S,c)
% find motif entries with duplicate nodes and correct by randomly swapping
% with another entry
N = length(S.V);
ix = find(S.I(:,c(1))==S.I(:,c(2)));
nstep = 100;
for n=1:nstep
if isempty(ix)
break;
end
for k=ix'
% select randomly to switch c(1) or (c2)
if rand<=0.5 % switch c(1)
col = c(1);
else % switch c(2)
col = c(2);
end
% randomly find entry for switching
n1 = ceil(N*rand);
for m=1:nstep
if ~ismember(S.I(k,col),S.I(n1,:)) && ~ismember(S.I(n1,col),S.I(k,:))
break;
end
n1 = ceil(N*rand);
end
if m==nstep
warning('Correctcolumns: maximum number of iterations reached (inner loop).');
end
% switch
S.I([k,n1],col) = S.I([n1,k],col);
end
ix = find(S.I(:,c(1))==S.I(:,c(2)));
end
if n==nstep
warning('Correctcolumns: maximum number of iterations reached (outer loop).');
end

function S = correctrows(S)
% find duplicate entries and correct by randomly swapping with another
% entry
N = length(S.V);
ndim = length(S.dim);
[I,Im,In] = unique(S.I,'rows');
nstep = 500;
for n=1:nstep
if length(I)==N
break;
end
% identify duplicate entries
for k=find(accumarray(In,ones(size(In)))>1)'
ix=find(In==k)';
% one entry can remain
for l=ix(2:end)
% pick random dimension
col = ceil(ndim*rand);
% randomly find entry for switching
n1 = ceil(N*rand);
for m=1:nstep
if ~ismember(S.I(l,col),S.I(n1,:)) && ~ismember(S.I(n1,col),S.I(l,:))
break;
end
n1 = ceil(N*rand);
end
if m==nstep
warning('Correctrows: maximum number of iterations reached (inner loop).');
end
% switch
S.I([l,n1],col) = S.I([n1,l],col);
end
end
[I,Im,In] = unique(S.I,'rows');
end
if n==nstep
warning('Correctcolumns: maximum number of iterations reached (outer loop).');
end

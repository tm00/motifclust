function [v,ev] = perronvector(varargin)
% PERRONVECTOR - Finds the Perron vector for a nonnegative symmetric matrix
% perronvector finds the Perron vector for a symmetric matrix s with
% nonnegative entries, i.e., the unique eigenvector with nonnegative
% entries corresponding to the largest eigenvalue. If s consists of
% disconnected components, the Perron vector for the largest component is
% computed. 
%
% USAGE:
%  [v,ev] = perronvector(s)
%
%    s  - symmetric matrix with nonnegative entries
%    ev - largest eigenvalue of s
%    v  - eigenvector corresponding to d with nonnegative entries
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

s = varargin{1};

% test symmetry
issym = size(s,1)==size(s,2) && isempty(find(s-s'));
if ~issym
error('Asymmetric input matrix is not allowed.')
end

% test for zero matrix
if nnz(s)==0
v = zeros(length(s),1);
ev = 0;
return;
end

opteigs.disp = 0;
opteigs.maxit = 300;
dim = length(s);

% we need nonzeros on the diagonal to get the disconnected components, see
% <http://blogs.mathworks.com/steve/2007/03/20/connected-component-labeling-part-3/>
if issparse(s)
s = s+speye(length(s));
else
s = s+eye(length(s));
end

[v,ev,flag] = eigs(s,2,'lm',opteigs);
if flag==0 && abs(ev(1,1)-ev(2,2))>1e-4
[ev, evi] = max(diag(ev));
v = v(:,evi);
else
% decompose s in disconnected components
[a,b,c,d] = dmperm(s);
if length(c)==2 % one component
if issparse(s) % use eigs for sparse matrix
[v,ev,flag] = eigs(s,1,'lm',opteigs);
if flag~=0
disp('eigs did not converge');
end
else % use eig for non-sparse matrix
[v,ev] = eig(s);
v = v(:,end);
ev = ev(end,end);
end
else % multiple components
sper = s(a,b); % block diagonally permuted s2
ev = 0;
% find largest component, if there are multiple components of the
% same size we need the one with highest eigenvalue
cc = c(2:end)-c(1:end-1);
blkidx = find(cc>1); %find(cc==max(cc));
if ~isempty(blkidx)
for ix = blkidx
%ix = blkidx(k);
% block matrix
blk = sper(c(ix):c(ix+1)-1,d(ix):d(ix+1)-1);
% 1st and last blocks don't correspond to separate component
% if not square
if size(blk,1)==size(blk,2)
if issparse(s) && length(blk)>1
[vtmp,etmp,flag] = eigs(blk,1,'lm',opteigs);
if flag~=0
disp('eigs did not converge');
end
else
[vtmp,etmp] = eig(full(blk));
vtmp = vtmp(:,end);
etmp = etmp(end,end);
end
% check if we increase emax
if etmp > ev
ev = etmp;
% embed block vector in whole space
v = zeros(dim,1);
v(c(ix):c(ix+1)-1) = vtmp;
% permute indices back to original order
[a2,ainv] = sort(a);
v = v(ainv);
end
end
end
else % s is diagonal matrix
[ev, imax] = max(diag(s));
v = zeros(dim,1);
v(imax) = 1;
end
end
end

% make sure v has nonnegative entries
[vmax, imax] = max(abs(v));
v = v./v(imax); 
% this one has all positive elements, upto numerical precision, hence:
v(v<10*eps) = 0;
v = sparse(v);
% return L2-normalized vector
v = v./norm(v);
% remember to subtract 1 from the eigenvalue
ev = ev-1;

function A = tensmat(varargin)
% TENSMAT - create matrix from tensor x vector inner product
% tensmat creates a matrix A from T and u by setting
%       d=1 => A(j,k) = sum_i u(i) T(i,j,k)
%       d=2 => A(i,k) = sum_j u(j) T(i,j,k)
%       d=3 => A(i,j) = sum_k u(k) T(i,j,k)
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
case 3 % sum over 1 dimension, return matrix
T = varargin{1};
u = varargin{2};
d = varargin{3};
% create new tensor T2(i,j,k) = u(i) T(i,j,k)
V = u(T.I(:,d)).*T.V;
T2.I = T.I(V>0,:);
T2.V = V(V>0);
T2.dim = T.dim;
% sum T2 over dimension d
dd = setdiff([1 2 3], d);
[I, m, n] = unique(T2.I(:,dd),'rows');
V = accumarray(n, T2.V);
I = I(V>0,:);
V = V(V>0);
% create A
A = sparse(I(:,1), I(:,2), V, T.dim(dd(1)), T.dim(dd(2)));
case 4 % sum over two dimensions, return vector
T = varargin{1};
u = varargin{2};
v = varargin{3};
d = varargin{4};
% create new tensor T2(i,j,k) = u(i) v(j) T(i,j,k)
V = u(T.I(:,d(1))).*v(T.I(:,d(2))).*T.V;
T2.I = T.I(V>0,:);
T2.V = V(V>0);
T2.dim = T.dim;
% sum T2 over dimension d(1) and d(2)
dd = setdiff([1 2 3], d);
[I, m, n] = unique(T2.I(:,dd));
V = accumarray(n, T2.V);
I = I(V>0);
V = V(V>0);
A = sparse(T.dim(dd),1);
A(I) = V;
case 5 % sum over three dimensions, return number
T = varargin{1};
u = varargin{2};
v = varargin{3};
w = varargin{4};
d = varargin{5};
% create new tensor T2(i,j,k) = u(i) v(j) w(k) T(i,j,k)
V = u(T.I(:,d(1))).*v(T.I(:,d(2))).*w(T.I(:,d(3))).*T.V;
A = sum(V);
end

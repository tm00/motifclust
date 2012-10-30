function [x,y,z,Sout,ovlp,nz] = tensor3clust(varargin)
% TENSOR3CLUST - Iterative clustering of a 3-dim adjacency tensor
% tensor3clust clusters a 3-dim tensor by computing its best rank 1
% approximation, interpreting them as a cluster weight, and using them to
% truncate the tensor for the next step.
%
% Tensor3clust is best called from wrapper functions which create the
% tensor for specific subgraph clustering problems
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

% proces input
switch nargin
    case 2
        T = varargin{1};
        status = varargin{2};
        p = 2;
    case 3
        T = varargin{1};
        status = varargin{2};
        p = varargin{3};
end

% initialize loop parameters
e = 2;
k = 0;
ev = [];
nz = [];
Sout = cell(1,1);
x = sparse(T.dim(1),1);
y = sparse(T.dim(2),1);
z = sparse(T.dim(3),1);
S.I = [];
S.V = [];
S.dim = T.dim;
ovlp = zeros(1,3);
while e>0 && size(T.I,1)>0
    k = k+1;
    % remove motifs in last cluster from adjacency tensor
    T = tensmin(T,S);
    if isempty(T.I)
        break;
    end
    % remove singletons
    %T = tensmin(T,tenssingl(T));
    % find maximizer for truncated adjacency tensor
    %[u,v,w,e] = maximizer3(T,status);
    %[u,v,w,e] = maximizer3knorm(T);
    [u,v,w,e] = maximizer3knorm(T,status,p);
    if isempty(find(u)) || isempty(find(v)) || isempty(find(w))
        break;
    end
    ev(k) = e;
    % get best uniform approximations
    if p==2
        [x(:,k), ovlp(k,1)] = unifapprox(u);
        [y(:,k), ovlp(k,2)] = unifapprox(v);
        [z(:,k), ovlp(k,3)] = unifapprox(w);
    else
        [x(:,k), ovlp(k,1)] = unifapproxknorm(u,p);
        [y(:,k), ovlp(k,2)] = unifapproxknorm(v,p);
        [z(:,k), ovlp(k,3)] = unifapproxknorm(w,p);
    end
    % create tensor
    S = tens3trunc(x(:,k),y(:,k),z(:,k),T);
    if size(S.I,1)==0 % rescue operation: take best 1-motif approximation
        warning('Used next best approximation for this cluster.')
        nz_resc = size(S.I,1);
        k_resc = 1;
        while nz_resc == 0
            k_resc = k_resc+1;
            [x(:,k), ovlp(k,1)] = unifapproxknorm(u,p,k_resc);
            [y(:,k), ovlp(k,2)] = unifapproxknorm(v,p,k_resc);
            [z(:,k), ovlp(k,3)] = unifapproxknorm(w,p,k_resc);
            [k_resc, nnz(x(:,k))]
            % create tensor
            S = tens3trunc(x(:,k),y(:,k),z(:,k),T);
            nz_resc = size(S.I,1);
        end
    end
    Sout{k} = S;
    nz(k) = size(S.I,1);
    % print some output
    disp(['CLUSTER ', num2str(k), ', nodes : (', ...
        num2str(length(unique(S.I(:,1)))), ', ', ...
        num2str(length(unique(S.I(:,2)))), ', ', ...
        num2str(length(unique(S.I(:,3)))), '), ev: ', num2str(e), ...
        ', nnz: (', num2str(size(S.I,1)), ', ', num2str(size(T.I,1)), ')']);
end
% sort output by eigenvalues
[ev, t] = sort(ev, 'descend');
x = x(:,t);
y = y(:,t);
z = z(:,t);
Sout = Sout(t);
ovlp = ovlp(t,:);
nz = nz(t);

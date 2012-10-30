function t = createtensor3(varargin)
% CREATETENSOR3 - create an adjacency tensor from a list of 3 adjacency matrices
% createtensor3 creates a sparse adjacency tensor t from a list of
% adjacency matrices A, B, C by setting
%
%          t(i,j,k) = A(i,j) B(j,k) C(k,i)
%
% (C can be omitted).
%
% t is a a structure with fields
%   I = non-zero indices
%   V = corresponding tensor values
%   dim = tensor dimensions
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de


switch nargin
    case 2
        A = varargin{1};
        B = varargin{2};
        % check if dimensions match
        if (size(A,2)~=size(B,1))
            error('Dimension mismatch.')
        end
        % for square matrices we want empty diagonals
        if size(A,1)==size(A,2)
            A = A - diag(diag(A));
        end
        if size(B,1)==size(B,2)
            B = B - diag(diag(B));
        end
        % start with empty tensor
        t.I = [];
        t.V = [];
        t.dim = [size(A,1) size(B,1) size(B,2)];
        % loop over nodes in C which participate in motif
        K = find(sum(A*B,1));
        for ix = 1:length(K)
            [I,J,val] = find(A.*(repmat(B(:,K(ix))',[size(A,1),1])));
            t.I = [t.I; I, J, repmat(K(ix),[length(I),1])];
            t.V = [t.V; val];
        end
        % we need i~=k
        NZ = find(t.I(:,1)~=t.I(:,3));
        t.I = t.I(NZ,:);
        t.V = t.V(NZ,:);
    case 3
        A = varargin{1};
        B = varargin{2};
        C = varargin{3};
        % check if dimensions match
        if (size(A,2)~=size(B,1) || size(B,2)~=size(C,1) || size(C,2)~=size(A,1))
            error('Dimension mismatch.')
        end
        if size(A,1)==size(A,2)
            A = A - diag(diag(A));
        end
        if size(B,1)==size(B,2)
            B = B - diag(diag(B));
        end
        if size(C,1)==size(C,2)
            C = C - diag(diag(C));
        end
        % start with empty tensor
        t.I = [];
        t.V = [];
        t.dim = [size(A,1) size(B,1) size(B,2)];
        % loop over nodes in C which participate in motif
        K = find(diag(C*A*B));
        for ix = 1:length(K)
            [I,J,val] = find(A.*(B(:,K(ix))*C(K(ix),:))');
            t.I = [t.I; I, J, repmat(K(ix),[length(I),1])];
            t.V = [t.V; val];
        end
    otherwise
        error('Wrong number of input arguments.')
end

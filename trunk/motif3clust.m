function varargout = motif3clust(varargin)
% MOTIF3CLUST - Iterative clustering of a 3-dim adjacency tensor
% motif3clust finds topological generalizations of the 3-node subgraphs
% by computing the best rank 1 approximation to the adjacency tensor,
% interpreting them as a cluster weight, and using them to truncate the
% adjacency tensor for the next step.
%
% USAGE:
%
%  The number of input/output arguments reflects symmetries in the motif:
%
%    (1) [x,y,z] = motif3clust(A1,A2,A3)
%        Most general case, no symmetries. This clusters the adjacency
%        tensor T(i,j,k) = A1(i,j)*A2(j,k)*A3(k,i)
%          A1 - sparse adjacency matrix for dimension (1,2)
%          A2 - sparse adjacency matrix for dimension (2,3)
%          A3 - sparse adjacency matrix for dimension (3,1)
%
%    (2) [x,y] = motif3clust(A1,A2)
%        Symmetry between dimensions (2,3). Only allowed if A2 is symmetric
%        and then means the same as motif3clust(A1,A2,A1') (but symmetry is
%        exploited in the algorithm).
%
%    (3) x = motif3clust(A1)
%        Symmetry between dimensions (1,2,3). Only allowed if A1 is symmetric
%        and then means the same as motif3clust(A1,A1,A1) (but symmetry is
%        exploited in the algorithm).
%
% OUTPUT:
%
%    (1) [x,y,z]
%       x  - matrix of cluster weights for dimension 1
%       y  - matrix of cluster weights for dimension 2
%       z  - matrix of cluster weights for dimension 3
%
%    (2) [x, y]
%       x  - matrix of cluster weights for dimension 1
%       y  - matrix of cluster weights for dimension 2 and 3
%
%    (3) x
%       x  - matrix of cluster weights for dimension 1, 2 and 3
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de/

% Process input and initialize output
switch nargin
    case 1 % 1 input matrix, should be symmetric
        A1 = varargin{1};
        issym1 = size(A1,1)==size(A1,2) && isempty(find(A1-A1'));
        if ~issym1
            error('1 input matrix requires symmetric matrix');
        end
        A1 = A1 - diag(diag(A1)); % we want empty diagonal
        % next 2 lines are faster than createtensor3(A1,A1,A1)
        T = createtensor3(triu(A1),triu(A1),triu(A1)');
        T = tenssym(T);
        nmTot = length(T.V)/6;
        status = 1;
    case 2 % 2 input matrices
        A1 = varargin{1};
        A2 = varargin{2};
        issym2 = size(A2,1)==size(A2,2) && isempty(find(A2-A2'));
        if ~issym2
            error('2 input matrices requires symmetric 2nd matrix');
        end
        if size(A1,1)==size(A1,2)
            A1 = A1 - diag(diag(A1)); % we want empty diagonal
        end
        if size(A2,1)==size(A2,2)
            A2 = A2 - diag(diag(A2));
        end
        T = createtensor3(A1,A2,A1');
        nmTot = length(T.V)/2;
        status = 2;
    case 3 % 3 input matrices
        A1 = varargin{1};
        A2 = varargin{2};
        A3 = varargin{3};
        if size(A1,1)==size(A1,2)
            A1 = A1 - diag(diag(A1)); % we want empty diagonal
        end
        if size(A2,1)==size(A2,2)
            A2 = A2 - diag(diag(A2));
        end
        if size(A3,1)==size(A3,2)
            A3 = A3 - diag(diag(A3));
        end
        T = createtensor3(A1,A2,A3);
        nmTot = length(T.V);
        status = 3;
    case 5 % same as 1 but matrix specified as edge lists
        A1 = sparse(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
        issym1 = size(A1,1)==size(A1,2) && isempty(find(A1-A1'));
        if ~issym1
            error('1 input matrix requires symmetric matrix');
        end
        A1 = A1 - diag(diag(A1)); % we want empty diagonal
        T = createtensor3(A1,A1,A1);
        nmTot = length(T.V)/6;
        status = 1;
    case 10 % same as 2 but matrices specified as edge lists
        A1 = sparse(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
        A2 = sparse(varargin{6},varargin{7},varargin{8},varargin{9},varargin{10});
        A3 = A1';
        issym2 = size(A2,1)==size(A2,2) && isempty(find(A2-A2'));
        if ~issym2
            error('2 input matrices requires symmetric 2nd matrix');
        end
        if size(A1,1)==size(A1,2)
            A1 = A1 - diag(diag(A1)); % we want empty diagonal
        end
        if size(A2,1)==size(A2,2)
            A2 = A2 - diag(diag(A2));
        end
        T = createtensor3(A1,A2,A1');
        nmTot = length(T.V)/2;
        status = 2;
    case 15 % same as 3 but matrices specified as edge lists
        A1 = sparse(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
        A2 = sparse(varargin{6},varargin{7},varargin{8},varargin{9},varargin{10});
        A3 = sparse(varargin{11},varargin{12},varargin{13},varargin{14},varargin{15});
        if size(A1,1)==size(A1,2)
            A1 = A1 - diag(diag(A1)); % we want empty diagonal
        end
        if size(A2,1)==size(A2,2)
            A2 = A2 - diag(diag(A2));
        end
        if size(A3,1)==size(A3,2)
            A3 = A3 - diag(diag(A3));
        end
        T = createtensor3(A1,A2,A3);
        nmTot = length(T.V);
        status = 3;
    otherwise
        error('Wrong number of input arguments.')
end

% check number of output variables
switch status
    case 1
        if nargout > 4
            error('Too many output arguments.')
        end
    case 2
        if nargout > 5
            error('Too many output arguments.')
        end
    case 3
        if nargout > 6
            error('Too many output arguments.')
        end
end

% cluster tensor
[x,y,z,Sout,ovlp,nz] = tensor3clust(T,status);

% output
switch status
    case 1
        % keep only modules with at least 2 motifs
        %sel = find(nz>6);
        sel = 1:length(nz);
        varargout{1} = x(:,sel);
        switch nargout
            case 2
                varargout{2} = Sout(sel);
            case 3
                varargout{2} = Sout(sel);
                varargout{3} = nmTot;
            case 4
                varargout{2} = Sout(sel);
                varargout{3} = nmTot;
                varargout{4} = ovlp(sel,1);
        end
    case 2
        % keep only modules with at least 2 motifs
        %sel = find(nz>2);
        sel = 1:length(nz);
        varargout{1} = x(:,sel);
        varargout{2} = y(:,sel);
        switch nargout
            case 3
                varargout{3} = Sout(sel);
            case 4
                varargout{3} = Sout(sel);
                varargout{4} = nmTot;
            case 5
                varargout{3} = Sout(sel);
                varargout{4} = nmTot;
                varargout{5} = ovlp(sel,1:2);
        end
    case 3
        % keep only modules with at least 2 motifs
        %sel = find(nz>1);
        sel = 1:length(nz);
        varargout{1} = x(:,sel);
        varargout{2} = y(:,sel);
        varargout{3} = z(:,sel);
        switch nargout
            case 4
                varargout{4} = Sout(sel);
            case 5
                varargout{4} = Sout(sel);
                varargout{5} = nmTot;
            case 6
                varargout{4} = Sout(sel);
                varargout{5} = nmTot;
                varargout{6} = ovlp(sel,:);
        end
end

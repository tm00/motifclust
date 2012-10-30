%% Network, tensor and motif representations
% This file lists some conventions regarding the representation of
% networks, tensors and motifs in the Network Motif Clustering Toolbox

%% Network representation
% Networks in Matlab are most easily represented as sparse adjacency
% matrices. Networks stored as lists of edges in text files can be imported
% using the <matlab:doc('importdata') importdata> function and converted to
% sparse matrices using the <matlab:doc('sparse') sparse> function. 
%
% A common scenario is to have an Nx2 matrix E of edges, where each row of
% E represents a source node and a target node. This is converted to an
% adjacency matrix using:
A = sparse(E(:,1),E(:,2),ones(length(E),1),max(E(:)),max(E(:)));
%%
% For a weighted network with edge weights in a third column of E, this
% becomes:
A = sparse(E(:,1),E(:,2),E(:,3),max(E(:)),max(E(:)));
%%
% By convention, all functions in the Network Motif Clustering Toolbox
% assume that an undirected network is given by a symmetric adjacency
% matrix and a directed network by an asymmetric adjacency matrix where
% A(i,j)>0 implies an edge from i to j. We also do not want self-edges to
% cause wrong motif statistics, so it is recommended to remove all self
% edges (if necessary, store self-edges for later use): 
Aself = find(diag(A));
A = A - diag(diag(A));
%% Tensor representation
% A tensor is a multi-dimensional array and we specifically need sparse
% tensors having few non-zero indices. In order to ensure compatibility of
% the Network Motif Clustering Toolbox with both <http://www.mathworks.com
% Matlab> and <http://www.gnu.org/software/octave/ Octave>, we use
% deliberately a very simple tensor representation. A tensor T is defined
% as a structure with 3 fields, I, V and dim. T.I is an Nx3 array
% containing the non-zero indices of T, S.V is an Nx1 vector of weights for
% the corresponding entries in S.I and S.dim are the total dimensions of
% the tensor. Thus a 10x10x10 tensor with non-zeros (1,1,1) and (1,2,3)
% is defined by
T.I = [1 1 1; 1 2 3];
T.V = [1; 1];
T.dim = [10 10 10];
%%
% A tensor class with a much richer feature set is availble in the
% <http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox/ MATLAB Tensor
% Toolbox>. Conversion to a sparse TensorToolbox tensor is trivial: 
S = sptensor(T.I,T.V,T.dim);
%% Motif representation
% The Network Motif Clustering Toolbox currently supports 3-node motifs.
% Throughout this documentation, a 3-node motif is represented as a string
% ABC. This representation means that if we cycle through the motif edges,
% the first edge comes from the network with adjacency matrix A, the second
% from B and the third from C. If for a directed network an edge is
% traversed in opposite direction, we need to use the matrix transpose. 
%%
% *Examples*
%
% * AAA = the triangle motif in an undirected or feedback loop motif in a
% directed network.
%
% * AAA' = the feedforward loop motif in a directed network
% 
% * ABA' = a mixed feedforward loop motif if B is directed or a symmetric
% "copointed" motif if B is undirected.
%
% * BAA' = a mixed feedforward loop motif if B is directed or a symmetric
% "copointing" motif if B is undirected.
%
% * etc.
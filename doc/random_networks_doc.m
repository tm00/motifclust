%% Random network ensembles
% This file describes the functions for generating random networks using
% the Network Motif Clustering Toolbox.

%% Randomized network with fixed degree distributions
% To generate a random network which preserves the in- and out-degree
% distribution of the input network by random edge swapping:
Ar = matrand(A);

%% Randomized network with fixed degree distributions and fixed number of motif instances
% The basic call for randomizing 3 networks A, B, C keeping the total
% number of instances of the network motif ABC (see 
% <network_representation_doc.html Network, tensor and motif
% representations>) constant, is:
[T,Ar] = tensor3rand(A,B,C,[1 2 3]);
%%
% The input arguments are a list of maximally 3 input networks to be
% randomized and a vector defining the motif. This vector cycles through the
% edges of a 3-node motif, the value indicates which network corresponds to
% each edge (according to the order of the input list), and a negative
% value indicates edge reversal, i.e. consistent with our
% <network_representation_doc.html 3-letter motif representations>. Some
% examples: 
%%
% The feedforward loop in a directed network A:
[T,Ar] = tensor3rand(A,[1 1 -1]);
%%
% Symmetric co-targeted motif with directed network A and undirected
% network B:
[T,Ar] = tensor3rand(A,B,[1 2 -1]);
%%
% Regarding the output arguments, T is the randomized tensor which can be
% used directly for input in the aggregation score calculation or
% clustering algorithm. Ar is a cell array of randomized networks, in the
% same order as the input arguments, so for the call
[T,Ar] = tensor3rand(A,B,C,[1 2 3]);
%%
% Ar{1} is the randomization of A, Ar{2} the randomization of B, etc.
%% 
% For clustering the feedforward loop in a single directed network, a
% slightly faster implementation with simpler user interface is available:
Ar = matrandffl(A);

%% Random network ensembles
% A random network ensemble is defined as a cell array of independently
% sampled random networks, for instance:
numNet = 1000;
Arand = cell(numNet,1);
for k=1:numNet
    Arand{k} = matrand(A);
end
%%
% For the tensor3rand function, we need the tensor to compute the
% aggregation statistic:
numNet = 1000;
Trand = cell(numNet,1);
for k=1:numNet
    Trand{k} = matrand(A,B,C,[1 2 3]);
end
%%
% It is wise to <matlab:doc('save') save> the ensemble for later use.

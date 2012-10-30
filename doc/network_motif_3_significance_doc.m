%% Network motif aggregation and enrichment significance
% This file describes how to use the network motif aggregation and
% enrichment functions of the Network Motif Clustering Toolbox for 3-node
% network motifs. 

%% Network motif aggregation significance
% For a given motif, compute its motif tensor T (using
% <matlab:doc('createtensor3') createtensor>) and generate and ensemble of
% randomized tensors Trand while keeping the degree distributions and total
% motif count of the corresponding input networks constant (see
% <random_networks_doc.html random network ensembles>). Then use:
[pval, Z, sreal, srand] = motif3aggregation(T, Trand);
%%
% where sreal and srand are the real and random aggregation statistics used
% to compute the P-value and Z-score.

%% Network motif enrichment significance
% Network motif enrichment can be easily computed for sets of motifs at a
% time. We distinguish between the number of 'edge colors' or integrated
% networks.
%%
% *One network*
[p,Z,nReal,nRand] = motif3enrich(A,Arand);
%% 
% *Two networks*
[p,Z,nReal,nRand] = motif3enrich(A,B,Arand,Brand);
%%
% *Three networks*
[p,Z,nReal,nRand] = motif3enrich(A,B,C,Arand,Brand,Crand);
%%
% *Remarks*
%
% * Arand is an ensemble of random networks (see <random_networks_doc.html
% random network ensembles>)
%
% * The function for two networks considers all motifs with at least one
% edge from each network, but not those with only one edge color that can
% be computed using the 'one network' variant. Idem for three networks
% which does not include the output that can be obtained with one or two
% networks. 
%
% * The number and type of possible motifs is determined by the (a)symmetry
% of the input matrices. It is assumed that if one matrix is symmetric, so
% are all the following, so care has to be given to the order of input
% arguments. The function does not given an error if for instance A is
% symmetric, but B is not, although the output will not be correct.
%
% * The output consists of vectors of p-value and Z-scores for each motif.
% To know which position corresponds to which motif, motif3enrich prints
% information on the screen when it is called. This information is also
% available in the help contents of motif3enrich and uses our
% <network_representation_doc.html standard motif representation
% conventions>. nReal is the number of  motifs in the real networks and
% nRand the number of motifs in each of the random networks. 
%
% * The function assumes that all real and random networks have empty
% diagonal. This is the case if random networks have been generated using
% the randomnet function. If necessary, remove diagonals using:
A = A - diag(diag(A));

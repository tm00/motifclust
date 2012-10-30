%% Network motif clustering
% This file describes how use the network motif clustering algorithm
% functions of the Network Motif Clustering Toolbox for 3-node network
% motifs. 

%% Calling the network motif clustering algorithm
% The basic call for clustering the network motif ABC (see
% <network_representation_doc.html Network, tensor and motif
% representations>) is: 
[x,y,z,S,o] = motif3clust(A,B,C);
nmTot = sum(diag(A*B*C));
%% 
% If the motif is symmetric for interchaning two nodes, for instance in the
% motif ABA' with B a symmetric matrix, then instead of using
% motif3clust(A,B,A'), use: 
[x,y,S,o] = motif3clust(A,B);
nmTot = sum(diag(A*B*A'))/2;
%% 
% For clustering the completely symmetric triangle motif in an undirected
% network, use:
[x,S,o] = motif3clust(A);
nmTot = sum(diag(A*A*A))/6;
%%
% We also listed each time the command to compute the total number of
% motifs (assuming all matrices have empty diagonal).

%% Output arguments
% *x, y, z - node sets*
%%
% For a motif ABC, x, y, z correspond to the motif nodes in the order such
% that the A-edge runs from x to y, the B-edge from y to z and the C-edge
% from z to x. x, y and z are matrices and each column corresponds to a
% cluster. Hence
X = find(x(:,k));
Y = find(y(:,k));
Z = find(z(:,k));
%%
% are the 3 sets of nodes which define cluster k.
%%
% *S - motif clusters*
%%
% S is a cell array of tensors where each tensor contains the motifs
% belonging to one cluster, see also <network_representation_doc.html
% Network, tensor and motif representations>. There is a certain redundancy
% in the output since  
X = unique(S{k}.I(:,1));
Y = unique(S{k}.I(:,2));
Z = unique(S{k}.I(:,3));
%%
% *o - overlap vectors*
%%
% In our papers we derived an upper bound for the distance between the
% motif cluster found using our algorithm and the true, unknown maximum.
% This upper bound goes to zero as the overlap between two vectors goes to
% one in each dimension (x, y, z). The output argument o contains these
% overlap values (|o(k,:)| are the overlap values for cluster k).

%% Post-processing functions
% The following functions compute useful quantities from the output of the
% motif clustering algorithm.
%%
% To compute the number of motifs and score of each cluster (the second
% command is only for the symmetric motif):
[nm,sc] = motif3cluststat(x,y,z,S);
[nm,sc] = motif3cluststat(x,y,S);

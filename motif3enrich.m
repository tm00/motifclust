function [p,Z,nReal,nRand] = motif3enrich(varargin)
% MOTIF3ENRICH - enrichment of a 3-node subgraph in random ensemble
% MOTIF3ENRICH computes the enrichment of a 3-node subgraph in an ensemble
% of random networks with the same expected degree distribution as the real
% network.
%
% USAGE: 
%
%  The number of input arguments reflects the number of integrated networks. 
%
%    (1) [p,Z,nReal,nRand] = motif3enrich(A,Arand)
%           - A sym: 1 motif (AAA)
%           - A asym: 2 motifs (AAA, AAA')
%
%    (2) [p,Z,nReal,nRand] = motif3enrich(A,B,Arand,Brand)
%           - A asym, B asym: 6 motifs (ABA, ABA', BAA', BAB, BAB', ABB')
%           - A asym, B sym: 4 motifs (ABA, ABA', BAA', ABB)
%           - A sym, B sym: 2 motifs (ABB, AAB)
%
%    (3) [p,Z,nReal,nRand] = motif3enrich(A,B,C,Arand,Brand,Crand)
%           - A asym, B asym, C asym: 8 motifs (ABC, ACB, ABC', ACB', BAC',
%             BCA', CAB', CBA')
%           - A asym, B asym, C sym: 4 motifs (ABC, ACB, ACB', CAB')
%           - A asym, B sym, C sym: 2 motifs (ABC, ACB)
%           - A sym, B sym, C sym: 2 motifs (ABC, ACB)
%
%  CONVENTIONS:
%     - It is assumed that all matrices (real and randomized) have empty
%       diagonal and the right symmetry properties.
%     - If a matrix is symmetric, it is assumed that all following matrices
%       are also symmetric.
%     - The random networks are given as a cell array of matrices with the
%       same dimension as the real network
%
%  OUTPUT PARAMETERS:
%     - p - enrichment p-value defined as the fraction of random networks
%       having at least the same number of subgraph instances as the real
%       network.
%     - Z - enrichment Z-score, (nReal-mean(nRand))/std(nRand).
%     - nReal - number of subgraphs in the real network.
%     - nRand - number of subgraphs in each of the random networks.
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
%  WRITTEN BY
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de

% Process input
switch nargin
case 2 % 1 input matrix
A = varargin{1};
Arand = varargin{2};
nE = length(Arand);
if issym(A)
disp('Symmetric matrix: AAA');
nReal = sum(diag(triu(A)*tril(A)*triu(A)));
nRand = zeros(nE,1);
for k=1:nE
nRand(k) = sum(diag(triu(Arand{k})*tril(Arand{k})*triu(Arand{k})));
end
else
disp('Asymmetric matrix: AAA, AAA''');
nReal = [sum(diag(A*A*A))/3, sum(diag(A*A*A'))];
nRand = zeros(nE,2);
for k=1:nE
nRand(k,1) = sum(diag(Arand{k}*Arand{k}*Arand{k}))/3;
nRand(k,2) = sum(diag(Arand{k}*Arand{k}*Arand{k}'));
end
end
case 4 % 2 input matrices
A = varargin{1};
B = varargin{2};
Arand = varargin{3};
Brand = varargin{4};
nE = length(Arand);
if ~issym(A) && ~issym(B)
disp('Two asymmetric matrices: ABA, ABA'', BAA'', BAB, BAB'', ABB''')
nReal = [sum(diag(A*B*A)), sum(diag(A*B*A')), sum(diag(B*A*A')), ...
sum(diag(B*A*B)), sum(diag(B*A*B')), sum(diag(A*B*B'))];
nRand = zeros(nE,6);
for k=1:nE
nRand(k,1) = sum(diag(Arand{k}*Brand{k}*Arand{k}));
nRand(k,2) = sum(diag(Arand{k}*Brand{k}*Arand{k}'));
nRand(k,3) = sum(diag(Brand{k}*Arand{k}*Arand{k}'));
nRand(k,4) = sum(diag(Brand{k}*Arand{k}*Brand{k}));
nRand(k,5) = sum(diag(Brand{k}*Arand{k}*Brand{k}'));
nRand(k,6) = sum(diag(Arand{k}*Brand{k}*Brand{k}'));
end
elseif ~issym(A) && issym(B)
disp('One asymmetric, one symmetric matrix: ABA, ABA'', BAA'', ABB')
nReal = [sum(diag(A*B*A)), sum(diag(A*B*A'))/2, ...
sum(diag(B*A*A'))/2, sum(diag(A*B*B'))];
nRand = zeros(nE,4);
for k=1:nE
nRand(k,1) = sum(diag(Arand{k}*Brand{k}*Arand{k}));
nRand(k,2) = sum(diag(Arand{k}*Brand{k}*Arand{k}'))/2;
nRand(k,3) = sum(diag(Brand{k}*Arand{k}*Arand{k}'))/2;
nRand(k,4) = sum(diag(Arand{k}*Brand{k}*Brand{k}'));
end
else
disp('Two symmetric matrices: AAB, ABB')
nReal = [sum(diag(A*A*B))/2, sum(diag(A*B*B))/2];
nRand = zeros(nE,2);
for k=1:nE
nRand(k,1) = sum(diag(Arand{k}*Arand{k}*Brand{k}))/2;
nRand(k,2) = sum(diag(Arand{k}*Brand{k}*Brand{k}))/2;
end
end
case 6 % 3 input matrices
A = varargin{1};
B = varargin{2};
C = varargin{3};
Arand = varargin{4};
Brand = varargin{5};
Crand = varargin{6};
nE = length(Arand);
if ~issym(A) && ~issym(B) && ~issym(C)
disp('Three asymmetric matrices: ABC, ACB, ABC'', ACB'', BAC'', BCA'', CAB'', CBA''')
nReal = [sum(diag(A*B*C)), sum(diag(A*C*B)), ...
sum(diag(A*B*C')), sum(diag(A*C*B')), ...
sum(diag(B*A*C')), sum(diag(B*C*A')), ...
sum(diag(C*A*B')), sum(diag(C*B*A'))];
nRand = zeros(nE,8);
for k=1:nE
nRand(k,1) = sum(diag(Arand{k}*Brand{k}*Crand{k}));
nRand(k,2) = sum(diag(Arand{k}*Crand{k}*Brand{k}));
nRand(k,3) = sum(diag(Arand{k}*Brand{k}*Crand{k}'));
nRand(k,4) = sum(diag(Arand{k}*Crand{k}*Brand{k}'));
nRand(k,5) = sum(diag(Brand{k}*Arand{k}*Crand{k}'));
nRand(k,6) = sum(diag(Brand{k}*Crand{k}*Arand{k}'));
nRand(k,7) = sum(diag(Crand{k}*Arand{k}*Brand{k}'));
nRand(k,8) = sum(diag(Crand{k}*Brand{k}*Arand{k}'));
end
elseif ~issym(A) && ~issym(B) && issym(C)
disp('Two asymmetric, one symmetric matrix: ABC, ACB, ACB'', CAB''')
nReal = [sum(diag(A*B*C)), sum(diag(A*C*B)), ...
sum(diag(A*C*B')), sum(diag(C*A*B'))];
nRand = zeros(nE,4);
for k=1:nE
nRand(k,1) = sum(diag(Arand{k}*Brand{k}*Crand{k}));
nRand(k,2) = sum(diag(Arand{k}*Crand{k}*Brand{k}));
nRand(k,3) = sum(diag(Arand{k}*Crand{k}*Brand{k}'));
nRand(k,4) = sum(diag(Crand{k}*Arand{k}*Brand{k}'));
end
elseif ~issym(A) && issym(B) && issym(C)
disp('One asymmetric, two symmetric matrices: ABC, ACB')
nReal = [sum(diag(A*B*C)), sum(diag(A*C*B))];
nRand = zeros(nE,2);
for k=1:nE
nRand(k,1) = sum(diag(Arand{k}*Brand{k}*Crand{k}));
nRand(k,2) = sum(diag(Arand{k}*Crand{k}*Brand{k}));
end
else
disp('Three symmetric matrices: ABC, ACB')
nReal = [sum(diag(A*B*C)), sum(diag(A*C*B))];
nRand = zeros(nE,2);
for k=1:nE
nRand(k,1) = sum(diag(Arand{k}*Brand{k}*Crand{k}));
nRand(k,2) = sum(diag(Arand{k}*Crand{k}*Brand{k}));
end
end
otherwise
error('Wrong number of input arguments.')
end


% p-value and Z-score
p = full(sum(nRand>=repmat(nReal, [nE,1]))./length(nRand));
% p = zeros(size(nReal));
% for k=1:length(nReal)
%     if nReal(k)<10
%         p(k) = full(sum(nRand(:,k)>=nReal(k)))/length(nRand(:,k));
%     else
%         [h,p(k)] = ztest(nReal(k), mean(nRand(:,k)), std(nRand(:,k)));
%     end
% end
Z = full((nReal-mean(nRand))./std(nRand));
nReal = full(nReal);
nRand = full(nRand);

function b = issym(A)
% ISSYM - check if input matrix is symmetric
b = size(A,1)==size(A,2) && isempty(find(A-A'));

function ar = matrandffl(a)
% MATRANDMOTIF - randomize network keeping degree distribution and total motif count constant
%
% This file is part of the Network Motif Clustering Toolbox
% Copyright 2011, Tom Michoel
% The full license terms can be found in Network_Motif_Clustering/LICENSE.txt
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de
%   http://omics.frias.uni-freiburg.de

a = a - diag(diag(a)); %remove diagonal elements
num = sum(diag(a*a*a')); % number of ffl
ar = matrand(a); % randomize keeping degree dist constant
numrand = sum(diag(ar*ar*ar')); % number of ffl

diff = 0.0;%0.01*num; % it's ok to be 1% off from real value
maxstep = 300;
step = 1;
while abs(numrand-num)>diff && step<maxstep
    disp(numrand-num)
    if numrand<num % randomly switch edges to "close" incomplete motifs
        % missing 1
        [I1,J1] = missingedges(ar,ar,ar');
        % missing 2
        [I2,J2] = missingedges(ar,ar',ar);
        % missing 3
        [I3,J3] = missingedges(ar',ar,ar);
        % sample numrand-num missing edges
        nn = length(I1) + length(I2) + length(I3);
        idm = ceil(rand(1,num-numrand)*nn);
        for id = idm(idm<=length(I1))
            i1 = I1(id);
            j1 = J1(id);
            if ar(i1,j1)==0
                K1 = find(ar(i1,:));
                j2 = ceil(rand*length(K1));
                K2 = find(ar(:,j1));
                i2 = ceil(rand*length(K2));
                % switch
                if j2>0 && i2>0
                    ar = edgeswap(ar,i1,K1(j2),K2(i2),j1);
                end
                % undo if we decreased number of motifs
                numr2 = sum(diag(ar*ar*ar'));
                if numr2<numrand
                    ar = edgeswap(ar,i1,j1,K2(i2),K1(j2));
                else
                    numrand = numr2;
                end
            end
        end
        for id = idm(idm>length(I1) & idm<=length(I1)+length(I2)) - length(I1)
            i1 = I2(id);
            j1 = J2(id);
            if ar(i1,j1)==0
                K1 = find(ar(i1,:));
                j2 = ceil(rand*length(K1));
                K2 = find(ar(:,j1));
                i2 = ceil(rand*length(K2));
                % switch
                if j2>0 && i2>0
                    ar = edgeswap(ar,i1,K1(j2),K2(i2),j1);
                end
                % undo if we decreased number of motifs
                numr2 = sum(diag(ar*ar*ar'));
                if numr2<numrand
                    ar = edgeswap(ar,i1,j1,K2(i2),K1(j2));
                else
                    numrand = numr2;
                end
            end
        end
        for id = idm(idm>length(I1)+length(I2)) - length(I1) - length(I2)
            i1 = J3(id);
            j1 = I3(id);
            if ar(i1,j1)==0
                K1 = find(ar(i1,:));
                j2 = ceil(rand*length(K1));
                K2 = find(ar(:,j1));
                i2 = ceil(rand*length(K2));
                % switch
                if j2>0 && i2>0
                    ar = edgeswap(ar,i1,K1(j2),K2(i2),j1);
                end
                % undo if we decreased number of motifs
                numr2 = sum(diag(ar*ar*ar'));
                if numr2<numrand
                    ar = edgeswap(ar,i1,j1,K2(i2),K1(j2));
                else
                    numrand = numr2;
                end
            end
        end
    else % randomly switch edges to "open" motifs
        T = createtensor3(ar,ar,ar');
        nn = length(T.V);
        idm = ceil(rand(numrand-num,1)*nn);
        for k=1:length(idm)
            id = idm(k);
            if ar(T.I(id,1),T.I(id,2))==1 && ar(T.I(id,2),T.I(id,3))==1 && ...
                    ar(T.I(id,1),T.I(id,3))==1
                % randomly pick edge
                ed = ceil(rand*3);
                switch ed
                    case 1 %  edge 1
                        [I,J] = find(ar);
                        jd = ceil(rand*length(I));
                        i1 = T.I(id,1);
                        j1 = T.I(id,2);
                        ar = edgeswap(ar,i1,j1,I(jd),J(jd));
                        % undo if we increased number of motifs
                        numr2 = sum(diag(ar*ar*ar'));
                        if numr2>numrand
                            ar = edgeswap(ar,i1,J(jd),I(jd),j1);
                        else
                            numrand = numr2;
                        end
                    case 2 %  edge 2
                        [I,J] = find(ar);
                        jd = ceil(rand*length(I));
                        i1 = T.I(id,2);
                        j1 = T.I(id,3);
                        ar = edgeswap(ar,i1,j1,I(jd),J(jd));
                        % undo if we increased number of motifs
                        numr2 = sum(diag(ar*ar*ar'));
                        if numr2>numrand
                            ar = edgeswap(ar,i1,J(jd),I(jd),j1);
                        else
                            numrand = numr2;
                        end
                    case 3 %  edge 3
                        [I,J] = find(ar);
                        jd = ceil(rand*length(I));
                        i1 = T.I(id,1);
                        j1 = T.I(id,3);
                        ar = edgeswap(ar,i1,j1,I(jd),J(jd));
                        % undo if we increased number of motifs
                        numr2 = sum(diag(ar*ar*ar'));
                        if numr2>numrand
                            ar = edgeswap(ar,i1,J(jd),I(jd),j1);
                        else
                            numrand = numr2;
                        end
                end
            end
        end
    end
    numrand = sum(diag(ar*ar*ar'));
    step = step+1;
end
if step==maxstep
    warning('tomic:tensor3rand', 'Tensor3rand: maximum number of steps reached, diff = %d', full(numrand-num));
end
%disp(step)

function [I,J] = missingedges(S1,S2,S3)
M = S2*S3;
[I1,J1] = find(S1.*M); % edges in motif
tmp = diag(sum(S1))*M*diag(sum(S1,2));
[J2,I2] = find(tmp-diag(diag(tmp)));
II = setdiff([I2,J2],[I1,J1],'rows');
I = II(:,1);
J = II(:,2);

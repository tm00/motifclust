function [T,Ar] = tensor3rand(varargin)
% TENSOR3RAND - generate random 3d tensor keeping network degree dist and total motif count constant

A = varargin(1:end-1);
mt = varargin{end};

Ar = cell(size(A));
for k=1:length(A)
    Ar{k} = matrand(A{k});
end
S = cell(3,1);
Sr = cell(3,1);
for k=1:length(S)
    S{k} = matselect(A{abs(mt(k))},mt(k));
    Sr{k} = matselect(Ar{abs(mt(k))},mt(k));
end

num = sum(diag(S{1}*S{2}*S{3})); % real number of motifs
numrand = sum(diag(Sr{1}*Sr{2}*Sr{3})); % random number of motifs

diff = 0.0; %0.01*num; % it's ok to be 1% off from real value
maxstep = 300;
step = 1;
while abs(numrand-num)>diff && step<maxstep
    disp(numrand-num)
    if numrand<num % randomly switch edges to "close" incomplete motifs
        % missing 1
        [I1,J1] = missingedges(Sr{1},Sr{2},Sr{3});
        % missing 2
        [I2,J2] = missingedges(Sr{2},Sr{3},Sr{1});
        % missing 3
        [I3,J3] = missingedges(Sr{3},Sr{1},Sr{2});
        % sample numrand-num missing edges
        nn = length(I1) + length(I2) + length(I3);
        nsample = min(num-numrand,500);
        idm = ceil(rand(1,nsample)*nn);
        %idm = ceil(rand(1,ceil(0.1*(num-numrand)))*nn);
        for id = idm(idm<=length(I1))
            % which direction
            if mt(1)>0
                i1 = I1(id);
                j1 = J1(id);
            else
                i1 = J1(id);
                j1 = I1(id);
            end
            if Ar{abs(mt(1))}(i1,j1)==0
                K1 = find(Ar{abs(mt(1))}(i1,:));
                j2 = ceil(rand*length(K1));
                K2 = find(Ar{abs(mt(1))}(:,j1));
                i2 = ceil(rand*length(K2));
                % switch
                if j2>0 && i2>0
                    Ar{abs(mt(1))} = edgeswap(Ar{abs(mt(1))},i1,K1(j2),K2(i2),j1);
                    Sr{1} = matselect(Ar{abs(mt(1))},mt(1));
                end
                % undo if we decreased number of motifs
                numr2 = sum(diag(Sr{1}*Sr{2}*Sr{3}));
                if numr2<numrand || numr2>num
                    Ar{abs(mt(1))} = edgeswap(Ar{abs(mt(1))},i1,j1,K2(i2),K1(j2));
                    Sr{1} = matselect(Ar{abs(mt(1))},mt(1));
                else
                    numrand = numr2;
                end
            end
        end
        for id = idm(idm>length(I1) & idm<=length(I1)+length(I2)) - length(I1)
            % which direction
            if mt(2)>0
                i1 = I2(id);
                j1 = J2(id);
            else
                i1 = J2(id);
                j1 = I2(id);
            end
            if Ar{abs(mt(2))}(i1,j1)==0
                K1 = find(Ar{abs(mt(2))}(i1,:));
                j2 = ceil(rand*length(K1));
                K2 = find(Ar{abs(mt(2))}(:,j1));
                i2 = ceil(rand*length(K2));
                % switch
                if j2>0 && i2>0
                    Ar{abs(mt(2))} = edgeswap(Ar{abs(mt(2))},i1,K1(j2),K2(i2),j1);
                    Sr{2} = matselect(Ar{abs(mt(2))},mt(2));
                end
                % undo if we decreased number of motifs
                numr2 = sum(diag(Sr{1}*Sr{2}*Sr{3}));
                if numr2<numrand  || numr2>num
                    Ar{abs(mt(2))} = edgeswap(Ar{abs(mt(2))},i1,j1,K2(i2),K1(j2));
                    Sr{2} = matselect(Ar{abs(mt(2))},mt(2));
                else
                    numrand = numr2;
                end
            end
        end
        for id = idm(idm>length(I1)+length(I2)) - length(I1) - length(I2)
            % which direction
            if mt(3)>0
                i1 = I3(id);
                j1 = J3(id);
            else
                i1 = J3(id);
                j1 = I3(id);
            end
            if Ar{abs(mt(3))}(i1,j1)==0
                K1 = find(Ar{abs(mt(3))}(i1,:));
                j2 = ceil(rand*length(K1));
                K2 = find(Ar{abs(mt(3))}(:,j1));
                i2 = ceil(rand*length(K2));
                % switch
                if j2>0 && i2>0
                    Ar{abs(mt(3))} = edgeswap(Ar{abs(mt(3))},i1,K1(j2),K2(i2),j1);
                    Sr{3} = matselect(Ar{abs(mt(3))},mt(3));
                end
                % undo if we decreased number of motifs
                numr2 = sum(diag(Sr{1}*Sr{2}*Sr{3}));
                if numr2<numrand || numr2>num
                    Ar{abs(mt(3))} = edgeswap(Ar{abs(mt(3))},i1,j1,K2(i2),K1(j2));
                    Sr{3} = matselect(Ar{abs(mt(3))},mt(3));
                else
                    numrand = numr2;
                end
            end
        end
    else % randomly switch edges to "open" motifs
        T = createtensor3(Sr{1},Sr{2},Sr{3});
        nn = length(T.V);
        nsample = numrand-num; %min(numrand-num,100);
        idm = ceil(rand(nsample,1)*nn);
        for k=1:length(idm)
            id = idm(k);
            if Sr{1}(T.I(id,1),T.I(id,2))==1 && Sr{2}(T.I(id,2),T.I(id,3))==1 && ...
                    Sr{3}(T.I(id,3),T.I(id,1))==1
                % randomly pick edge
                ed = ceil(rand*3);
                switch ed
                    case 1 %  edge 1
                        [I,J] = find(Ar{abs(mt(1))});
                        jd = ceil(rand*length(I));
                        if mt(1)>0
                            i1 = T.I(id,1);
                            j1 = T.I(id,2);
                        else
                            i1 = T.I(id,2);
                            j1 = T.I(id,1);
                        end
                        Ar{abs(mt(1))} = edgeswap(Ar{abs(mt(1))},i1,j1,I(jd),J(jd));
                        Sr{1} = matselect(Ar{abs(mt(1))},mt(1));
                        % undo if we increased number of motifs
                        numr2 = sum(diag(Sr{1}*Sr{2}*Sr{3}));
                        if numr2>numrand || numr2<num
                            Ar{abs(mt(1))} = edgeswap(Ar{abs(mt(1))},i1,J(jd),I(jd),j1);
                            Sr{1} = matselect(Ar{abs(mt(1))},mt(1));
                        else
                            numrand = numr2;
                        end
                    case 2 %  edge 2
                        [I,J] = find(Ar{abs(mt(2))});
                        jd = ceil(rand*length(I));
                        if mt(2)>0
                            i1 = T.I(id,2);
                            j1 = T.I(id,3);
                        else
                            i1 = T.I(id,3);
                            j1 = T.I(id,2);
                        end
                        Ar{abs(mt(2))} = edgeswap(Ar{abs(mt(2))},i1,j1,I(jd),J(jd));
                        Sr{2} = matselect(Ar{abs(mt(2))},mt(2));
                        % undo if we increased number of motifs
                        numr2 = sum(diag(Sr{1}*Sr{2}*Sr{3}));
                        if numr2>numrand || numr2<num
                            Ar{abs(mt(2))} = edgeswap(Ar{abs(mt(2))},i1,J(jd),I(jd),j1);
                            Sr{2} = matselect(Ar{abs(mt(2))},mt(2));
                        else
                            numrand = numr2;
                        end
                    case 3 %  edge 3
                        [I,J] = find(Ar{abs(mt(3))});
                        jd = ceil(rand*length(I));
                        if mt(3)>0
                            i1 = T.I(id,3);
                            j1 = T.I(id,1);
                        else
                            i1 = T.I(id,1);
                            j1 = T.I(id,3);
                        end
                        Ar{abs(mt(3))} = edgeswap(Ar{abs(mt(3))},i1,j1,I(jd),J(jd));
                        Sr{3} = matselect(Ar{abs(mt(3))},mt(3));
                        % undo if we increased number of motifs
                        numr2 = sum(diag(Sr{1}*Sr{2}*Sr{3}));
                        if numr2>numrand || numr2<num
                            Ar{abs(mt(3))} = edgeswap(Ar{abs(mt(3))},i1,J(jd),I(jd),j1);
                            Sr{3} = matselect(Ar{abs(mt(3))},mt(3));
                        else
                            numrand = numr2;
                        end
                end
            end
        end
    end
    Sr{1} = matselect(Ar{abs(mt(1))},mt(1));
    Sr{2} = matselect(Ar{abs(mt(2))},mt(2));
    Sr{3} = matselect(Ar{abs(mt(3))},mt(3));
    numrand = sum(diag(Sr{1}*Sr{2}*Sr{3}));
    step = step+1;
end
T = createtensor3(Sr{1},Sr{2},Sr{3});
if step==maxstep
    warning('tomic:tensor3rand', 'Tensor3rand: maximum number of steps reached, diff = %d', full(numrand-num));
end
%disp(step)

function s = matselect(A,t)
% A cell array, t index, sign t for transpose
if t>0
    s = A;
elseif t<0
    s = A';
else
    error('index must be nonzero')
end

function [I,J] = missingedges(S1,S2,S3)
M = S2*S3;
[I1,J1] = find(S1.*M); % edges in motif
tmp = diag(sum(S1))*M*diag(sum(S1,2));
[J2,I2] = find(tmp-diag(diag(tmp)));
II = setdiff([I2,J2],[I1,J1],'rows');
I = II(:,1);
J = II(:,2);

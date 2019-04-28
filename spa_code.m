function [blocks,D,idx]=spa_code(data,sigma,D)
addpath('');    % the file path of KSVD_Matlab_Toolbox should be added
reduceDC = 0;
[NN1,NN2] = size(data);
if (sigma > 5)
    numIterOfKsvd = 5;   % iteration number 5
else
    numIterOfKsvd = 5;
end
C = 1.15;
maxBlocksToConsider = 280000;   % 260000
slidingDis = 1;
bb = 8;
maxNumBlocksToTrainOn = 65000;
displayFlag = 1;
if (sigma <= 5)
    numIterOfKsvd = 5;
end
% first, train a dictionary on blocks from the noisy image
if(prod([NN1,NN2]-bb+1)> maxNumBlocksToTrainOn)
    randPermutation =  randperm(prod([NN1,NN2]-bb+1));
    selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);
    blkMatrix = zeros(bb^2,maxNumBlocksToTrainOn);
    for i = 1:maxNumBlocksToTrainOn
        [row,col] = ind2sub(size(data)-bb+1,selectedBlocks(i));
        currBlock = data(row:row+bb-1,col:col+bb-1);
        blkMatrix(:,i) = currBlock(:);
    end
else
    blkMatrix = im2col(data,[bb,bb],'sliding');
end

param.K = 256;
param.numIteration = numIterOfKsvd ;
param.errorFlag = 1; % decompose signals until a certain error is reached. do not use fix number of coefficients.
param.errorGoal = sigma*C;
param.preserveDCAtom = 0;

Pn=ceil(sqrt(param.K));
DCT=zeros(bb,Pn);
for k=0:1:Pn-1,
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);
if D==0
    param.initialDictionary = DCT(:,1:param.K );
else
    param.initialDictionary=D;
end
param.InitializationMethod =  'GivenMatrix';

if (reduceDC)
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end
param.displayProgress = displayFlag;
%[Dictionary,output] = KSVD(blkMatrix,param);   % if there is no
%dictionary, this line should not be skipped.
load Dictionary    % load a trained dictionary, if there is no dictionary, this line should be skipped

if (displayFlag)
    disp('finished Trainning dictionary');
end
% denoise the image using the resulted dictionary
errT = sigma*C;
% while (prod(floor((size(data)-bb)/slidingDis)+1)>maxBlocksToConsider)
%     slidingDis = slidingDis+1;
% end
[blocks,idx] = my_im2col(data,[bb,bb],slidingDis);
blk=parallsparsecoding(blocks,Dictionary,errT);
blocks=blk;
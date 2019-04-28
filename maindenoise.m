function datadenoised=maindenoise(data,B,sigma)

% regularization parameters, they should be tuned to obtain the best performance
lambda1=1;   % regularization parameter of non-local
lambda2=1;   % regularization parameter of low rank
mu1=1;    % parameter of auxiliry viarible
mu2=1;
iteration_num=15;
X=data;
D=0;   % the dictionary is initialized 0

for itr=1:iteration_num
    
    % optimization with respect of alpha and D using sparse coding
    [cA cH cV cD]=dwt2(X,'db2');
    Nsig=median(abs(cD(:)))/0.6745;
    [blocks,D,idx]=spa_code(X,Nsig,D);

    % optimization with respect of U, if non-local regularizer is not used,
    % it can be skipped
    R  =  1:size(B,1);
    C  =  1:size(B,2);
    V  =  1;
    mI =  sparse(R,C,V,size(B,1),size(B,2));
    U=(lambda1*B'*B+mu1*mI)\(mu1*X);

    % optimization with respect of V
    [S1 S2 S3]=svd(X,'econ');
    V=S1*soft(S2,lambda2/(2*mu2))*S3';
% 
%     % optimization with respect of X
    bb=8;
    count = 1;
    Weight = zeros(size(data));
    IMout = zeros(size(data));
    [rows,cols] = ind2sub(size(data)-bb+1,idx);
    for i  = 1:length(cols)
        col = cols(i); row = rows(i);        
        block =reshape(blocks(:,count),[bb,bb]);
        IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;
        Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);
        count = count+1;
    end;
    X = ((30/sigma)*data+IMout+mu1*U+mu2*V)./(30/sigma+Weight+mu1+mu2);

end
datadenoised=X;
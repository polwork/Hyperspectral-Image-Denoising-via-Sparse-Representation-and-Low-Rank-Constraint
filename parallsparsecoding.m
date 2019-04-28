function blocks2=parallsparsecoding(blocks,Dictionary,errT)
blocks2=zeros(size(blocks));
parfor jj = 1:size(blocks,2)   % 30000
    
    if (0)
        vecOfMeans = mean(blocks(:,jj));
        blocks2(:,jj) = blocks(:,jj) - repmat(vecOfMeans,size(blocks,1),1);
    end
    Coefs = OMPerr(Dictionary,blocks(:,jj),errT);
    if (0)
        %blocks2(:,jj)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
    else
        blocks2(:,jj)= Dictionary*Coefs ;
    end
end

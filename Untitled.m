% this code is the main code, which is to denoise the hyperspectral imagery
% on the assumption of additive Gaussian noise. It denoise in the framework
% of sparsity and low rank. Details of the method is presented in the paper
% published on IGARSS2013.
% Author: Jingxiang Yang
% Date: March 27 2013

clear all
e='';    % the file type of the HSI image, e.g., '.bmp','.jpg'
s='';    % the file path of the HSI image

numband=224;        %number of bands
row=256;         % number of rows of the HSI
col=256;         % number of the columns of the HSI
sigma=0.10;          %standard of variance of Gaussian noise


originaldata=zeros(row,col,numband);
noisedata=zeros(row,col,numband);
data=zeros(row*col,numband);    %size of HSI must be 75*75
fildata=zeros(row,col,numband);
PSNRIn=zeros(1,numband);
PSNROut=zeros(1,numband);
SSIMIn=zeros(1,numband);
FSIMIn=zeros(1,numband);
SSIMOut=zeros(1,numband);
FSIMOut=zeros(1,numband);
n=1;

for num=1:1:numband    % read HSI image
    m=int2str(num);
    name=strcat(s,m,e);
    I=double(imread(name));
    originaldata(:,:,num)=I;
    II=I+sigma*randn(size(I));  % Gaussian noise added to HSI
    noisedata(:,:,num)=II;
    data(:,n)=II(:);
    n=n+1;
end

save originaldata originaldata
save noisedata noisedata
%B=Compute_NLM_Matrix(originaldata(:,:,10),5);   % calculate non-local similarity matrix in spatial dimension, if non-local regularizer is not used, B can be set 0
B=0;

datadenoised=maindenoise(data,B,sigma);    % begin to denoise

fildata=reshape(datadenoised,row,col,numband);

for num=1:numband    % compute the assessing indices
    PSNRIn(1,num) = 20*log10(255/sqrt(mean(mean((noisedata(:,:,num)-originaldata(:,:,num)).^2))));
    PSNROut(1,num) = 20*log10(255/sqrt(mean(mean((fildata(:,:,num)-originaldata(:,:,num)).^2))));
    SSIMIn(1,num)    =    cal_ssim(noisedata(:,:,num),originaldata(:,:,num),0,0);
    [fism1, FSIMc] =    FeatureSIM(noisedata(:,:,num),originaldata(:,:,num));
    FSIMIn(1,num)=fism1;
    SSIMOut(1,num)    =    cal_ssim(fildata(:,:,num),originaldata(:,:,num),0,0);
    [fism2, FSIMc] =    FeatureSIM(fildata(:,:,num),originaldata(:,:,num));
    FSIMOut(1,num)=fism2;
end
MPSNRIn=mean(PSNRIn);
MPSNROut=mean(PSNROut);
MSSIMIn=mean(SSIMIn);
MSSIMOut=mean(SSIMOut);
MFSIMIn=mean(FSIMIn);
MFSIMOut=mean(FSIMOut);
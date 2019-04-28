function  [mW] =  Compute_NLM_Matrix( im, ws )
%----------------------------
% Only for grayscale image
% Apr. 13, 2010
%----------------------------

S       =  12;
f       =  ws;
t       =  floor(f/2);
nv      =  10;  %par.nblk;
hp      =  65;

e_im    =  padarray( im, [t t], 'symmetric' );
[h w]   =  size( im );
nt      =  (nv)*h*w;
R       =  zeros(nt,1);
C       =  zeros(nt,1);
V       =  zeros(nt,1);

L       =  h*w;
X       =  zeros(f*f, L, 'single');

% For the Y component
k    =  0;
for i  = 1:f
    for j  = 1:f
        k    =  k+1;
        blk  =  e_im(i:end-f+i,j:end-f+j);
        X(k,:) =  blk(:)';
    end
end

% Index image
I    =   reshape((1:L), h, w);
X    =   X'; 
f2   =   f^2;

cnt     =  1;
for  row  =  1 : h
    for  col  =  1 : w
        
        off_cen  =  (col-1)*h + row;
        
        rmin    =   max( row-S, 1 );
        rmax    =   min( row+S, h );
        cmin    =   max( col-S, 1 );
        cmax    =   min( col+S, w );
         
        idx     =   I(rmin:rmax, cmin:cmax);
        idx     =   idx(:);
        B       =   X(idx, :);        
        v       =   X(off_cen, :);
        
        
        dis     =   (B(:,1) - v(1)).^2;
        for k = 2:f2
            dis   =  dis + (B(:,k) - v(k)).^2;
        end
        dis   =  dis./f2;
        [val,ind]   =  sort(dis);        
        dis(ind(1))  =  dis(ind(2));        
        wei         =  exp( -dis(ind(1:nv))./hp );
        
        R(cnt:cnt+nv-1)     =  off_cen;
        C(cnt:cnt+nv-1)     =  idx( ind(1:nv) );
        V(cnt:cnt+nv-1)     =  wei./(sum(wei)+eps);
        cnt                 =  cnt + nv;
        
    end
end
R    =  R(1:cnt-1);
C    =  C(1:cnt-1);
V    =  V(1:cnt-1);
W1   =  sparse(R, C, V, h*w, h*w);

R    =  zeros(h*w,1);
C    =  zeros(h*w,1);
V    =  zeros(h*w,1);

R(1:end)  =  1:h*w;
C(1:end)  =  1:h*w;
V(1:end)  =  1;
mI        =  sparse(R,C,V,h*w,h*w);
mW        =  mI - W1;





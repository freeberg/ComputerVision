function [ncc,outside_image] = compute_ncc(d,im_right,P_right,im_left,segm_left,P_left,patchsize,sc)

%Change coordinates so that P_right = K[I 0];
[K,R] = rq(P_right);
T = [R(:,1:3)' -R(:,1:3)'*R(:,4); 0 0 0 1];
P_right = P_right*T;
P_left = P_left*T;

%rescale images
im_right = imresize(im_right,sc);
im_left = imresize(im_left,sc);
segm_left = imresize(segm_left,sc);
[K,R] = rq(P_left);
K(1:2,:) = K(1:2,:)*sc;
P_left = K*R;
[K,R] = rq(P_right);
K(1:2,:) = K(1:2,:)*sc;
P_right = K*R;

disp('Computing cross correlations:');
ncc = zeros(size(im_right,1),size(im_right,2),length(d));
outside_image = false(size(ncc));

%compute ncc-mean and ncc-norm of right image
meanpatch = ones(2*patchsize+1)./sum(sum(ones(2*patchsize+1)))/3;
patch = ones(2*patchsize+1);
Rright = double(im_right(:,:,1));
Gright = double(im_right(:,:,2));
Bright = double(im_right(:,:,3));
mean_right = imfilter(Rright,meanpatch,'same')+imfilter(Gright,meanpatch,'same')+imfilter(Bright,meanpatch,'same');

term1R = imfilter(Rright.^2,patch,'same');
term1G = imfilter(Gright.^2,patch,'same');
term1B = imfilter(Bright.^2,patch,'same');

term2R = mean_right.*imfilter(Rright,patch,'same');
term2G = mean_right.*imfilter(Gright,patch,'same');
term2B = mean_right.*imfilter(Bright,patch,'same');


term4 = sum(patch(:))*3*mean_right.^2;
norm_right = sqrt(term1R+term1G+term1B-2*(term2R+term2G+term2B)+term4);

for i = 1:length(d);
    %transform left image by current disp
    fprintf('Depth %d of %d\n',i,length(d));
    K1 = rq(P_right);
    [K2,R2] = rq(P_left);
    pi = [0 0 1 -d(i)];
    H = inv(K2*(R2(:,1:3)*pi(4)-R2(:,4)*pi(1:3))*inv(K1));
    tform = maketform('projective',H');
    imtr = imtransform(im_left,tform,'xdata',[1 size(im_right,2)],'ydata',[1 size(im_right,1)],'size',size(im_right));
    segmtr = imtransform(segm_left,tform,'xdata',[1 size(im_right,2)],'ydata',[1 size(im_right,1)],'size',size(im_right(:,:,1)));
    org_bnd_im = ones(size(im_left(:,:,1)));
    %To remove edge effects
    org_bnd_im(:,1:patchsize) = 0;
    org_bnd_im(:,end-patchsize:end) = 0;
    
    bnd_im = imtransform(org_bnd_im,tform,'xdata',[1 size(im_right,2)],'ydata',[1 size(im_right,1)],'size',size(im_right(:,:,1)),'fill',0);
    if 1
        figure(1);
        imagesc(uint8(abs(double(im_right)-double(imtr))))
    end

    Rtr = double(imtr(:,:,1));
    Gtr = double(imtr(:,:,2));
    Btr = double(imtr(:,:,3));
    
    mean_tr = imfilter(Rtr,meanpatch,'same')+imfilter(Gtr,meanpatch,'same')+imfilter(Btr,meanpatch,'same');
    %Compute ncc-mean and ncc- norm of transformed-left image
    
    term1R = imfilter(Rtr.^2,patch,'same');
    term1G = imfilter(Gtr.^2,patch,'same');
    term1B = imfilter(Btr.^2,patch,'same');
    
    term2R = mean_tr.*imfilter(Rtr,patch,'same');
    term2G = mean_tr.*imfilter(Gtr,patch,'same');
    term2B = mean_tr.*imfilter(Btr,patch,'same');

    
    term4 = sum(patch(:))*3*mean_tr.^2;
    norm_tr = sqrt(term1R+term1G+term1B-2*(term2R+term2G+term2B)+term4);

    %Compute ncc
    
    term1R = imfilter(Rright.*Rtr,patch,'same');
    term1G = imfilter(Gright.*Gtr,patch,'same');
    term1B = imfilter(Bright.*Btr,patch,'same');
    
    term2R = mean_right.*imfilter(Rtr,patch,'same');
    term2G = mean_right.*imfilter(Gtr,patch,'same');
    term2B = mean_right.*imfilter(Btr,patch,'same');
    
    term3R = mean_tr.*imfilter(Rright,patch,'same');
    term3G = mean_tr.*imfilter(Gright,patch,'same');
    term3B = mean_tr.*imfilter(Bright,patch,'same');
    term4 = sum(patch(:))*3*mean_tr.*mean_right;
    
    ncci = term1R+term1G+term1B-(term2R+term2G+term2B)-(term3R+term3G+term3B)+term4;
    ncci = ncci./norm_right./norm_tr;
    
    ncci(~isfinite(ncci)) = 0;
    ncci(~segmtr) = -1;
    ncci(~(bnd_im>=1-1e-8)) = 0;

    outside_image(:,:,i) = ~(bnd_im>=1-1e-8);
    ncc(:,:,i) = real(ncci);
    if 1
        figure(2);
        imagesc(real(ncci));
        caxis([-1 1]);
        colorbar
        %pause;
    end
end


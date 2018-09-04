clear
x = {};
im1 = imread('a.jpg');
im2 = imread('b.jpg');

points1 = detectSURFFeatures(rgb2gray(im1));
points2 = detectSURFFeatures(rgb2gray(im2));
imshow(im1)
hold on
plot(points1.selectStrongest (200));
%%
[features1 , validPoints1] = extractFeatures(rgb2gray(im1),points1 );
[features2 , validPoints2] = extractFeatures(rgb2gray(im2),points2 );
matches = matchFeatures(features1 , features2 , 'Unique', true);
x{1} = validPoints1(matches(:, 1)). Location';
x{2} = validPoints2(matches(:, 2)). Location';
[n, kx] = size(x{1});

x{1} = [x{1}; ones(1,kx)];
x{2} = [x{2}; ones(1,kx)];
%%
perm = randperm(size(matches ,1));
figure;
imagesc ([im1 im2]);
hold on;
plot([x{1}(1,perm (1:10));  x{2}(1,perm (1:10))+ size(im1 ,2)], ...
[x{1}(2,perm (1:10));  x{2}(2,perm (1:10))] ,'-');
hold  off;
%%

[n, kx] = size(x{1})
bestInl = 0;
bestH = []
bestx1H = []

for i=1:10
  perm = randperm(kx, 8);
  M=[];
  for col=1:8
    for row=1:n
        M=[M; zeros(1,3*(row-1)), x{2}(:,perm(col))', zeros(1,3*(3-row)),zeros(1,col-1),x{1}(row, perm(col)),zeros(1,(8-col))];
    end
  end
  [U,S,V] = svd(M);
  solu = V(1:9, end);
  H = [solu(1:3)'; solu(4:6)'; solu(7:9)'];
  
  x1H = H * x{2};
  x1H = pflat(x1H);
  inliers = (sqrt((x1H(1,:) - x{1}(1,:)).^2 + (x1H(2,:) - x{1}(2,:)).^2) <= 5);
  amtInl = sum(inliers(:) == 1);
  if (amtInl > bestInl)
       bestInl = amtInl;
       bestInliers = inliers;
       bestH = H;
       bestx1H = x1H;
   end
end

%%
bestH = double(bestH);
tform = maketform('projective', bestH');
transfbounds = findbounds(tform ,[1 1; size(im2,2) size(im2 ,1)]);
xdata = [min([ transfbounds(:,1); 1]) max([ transfbounds(:,1);  size(im2 ,2)])];
ydata = [min([ transfbounds(:,2); 1]) max([ transfbounds(:,2);  size(im2 ,1)])];
[newA] = imtransform(im2,tform ,'xdata',xdata ,'ydata',ydata);
imagesc(newA);
%%
tform2 = maketform('projective',eye (3));
transfbounds = findbounds(tform ,[1 1; size(im1,2) size(im1 ,1)]);
xdata = [min([ transfbounds(:,1); 1]) max([ transfbounds(:,1);  size(im1 ,2)])];
ydata = [min([ transfbounds(:,2); 1]) max([ transfbounds(:,2);  size(im1 ,1)])];
[newB] = imtransform(im1,tform2 ,'xdata',xdata ,'ydata',ydata ,'size',size(newA));

newAB = newB;
newAB(newB < newA) = newA(newB < newA);
imagesc(newAB);
%%

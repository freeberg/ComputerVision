%Exercise 4
clear
cube1 = imread('cube1.JPG');
cube2 = imread('cube2.JPG');

points1 = detectSURFFeatures(rgb2gray(cube1));
points2 = detectSURFFeatures(rgb2gray(cube2));

figure
imshow(cube1);
hold on
plot(points1.selectStrongest (200));
figure
imshow(cube2);
hold on
plot(points2.selectStrongest(200));

[features1, validPoints1] = extractFeatures(rgb2gray(cube1),points1 );
[features2, validPoints2] = extractFeatures(rgb2gray(cube2),points2 );

matches = matchFeatures(features1 , features2 , 'Unique', true);

x1 = validPoints1(matches(:, 1)).Location';
x2 = validPoints2(matches(:, 2)).Location';

perm = randperm(size(matches,1));
figure;
imagesc ([cube1 cube2]);
hold on;
plot([x1(1,perm(1:10));  x2(1,perm(1:10))+ size(cube1,2)], ...
     [x1(2,perm(1:10));  x2(2,perm(1:10))],'-');
hold  off;
%%
%Exercise 5
load('cameraP1.mat');
load('cameraP2.mat');
[n, kx] = size(x1)

x1h = [x1; ones(1,kx)];
x2h = [x2; ones(1,kx)];

DLT={};
eqnbr = 0;
for i=1:kx
    eqnbr = eqnbr+1;
    DLT{eqnbr}=[P1, x1h(:,i), zeros(3, 1); P2, zeros(3, 1), x2h(:,i)];
end

X3D = []
for j=1:kx
    [U,S,V] = svd(DLT{j});
    solu = V(:,end);
    X3D = [X3D, solu(1:4,1)];
end

xproj1 = pflat(P1*X3D);
xproj2 = pflat(P2*X3D);

figure
imagesc(cube1);
hold on
plot(xproj1(1,:), xproj1(2,:), '*')
plot(x1h(1,:), x1h(2,:), 'o')

figure
imagesc(cube2);
hold on
plot(xproj2(1,:), xproj2(2,:), '*')
plot(x2h(1,:), x2h(2,:), 'o')

%%
%Compare when normalized with inv K
[K1 Q1] = rq(P1);
[K2 Q2] = rq(P2);

xproj1k = inv(K1) * xproj1;
xproj2k = inv(K2) * xproj2;
x1hk = inv(K1) * x1h;
x2hk = inv(K2) * x2h;


figure
plot(xproj1k(1,:), xproj1k(2,:), '*')
hold on
plot(x1hk(1,:), x1hk(2,:), 'o')

figure
plot(xproj2k(1,:), xproj2k(2,:), '*')
hold on
plot(x2hk(1,:), x2hk(2,:), 'o')

%%
%Compute the errors
good_points = (sqrt(sum((x1-xproj1(1:2 ,:)).^2))  < 3 & ...
               sqrt(sum((x2-xproj2(1:2 ,:)).^2))  < 3);
           
Xgp = X3D(:,good_points);
Xgp = pflat(Xgp);
x1gp = x1h(:,good_points);
x2gp = x2h(:,good_points);

xproj1gp = pflat(P1*Xgp);
xproj2gp = pflat(P2*Xgp);

figure
imagesc(cube1);
hold on
plot(xproj1gp(1,:), xproj1gp(2,:), '*')
plot(x1gp(1,:), x1gp(2,:), 'o')

figure
imagesc(cube2);
hold on
plot(xproj2gp(1,:), xproj2gp(2,:), '*')
plot(x2gp(1,:), x2gp(2,:), 'o')

%3d plot
load('compEx3data.mat');
figure
plot3(Xgp(1,:), Xgp(2,:), Xgp(3,:), '.')
hold on
plot3([Xmodel(1,startind);  Xmodel(1,endind)],...
      [Xmodel(2,startind );  Xmodel(2,endind )],...
      [Xmodel(3,startind );  Xmodel(3,endind)],'b-');
cams = {1,2};
cams{1} = P1;
cams{2} = P2;
plotcams(cams)




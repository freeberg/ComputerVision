clear
load EssenMatrix.mat
load FundalMatrix.mat
load compEx1data.mat
load compEx3data.mat

W = [0 -1 0; 1 0 0; 0 0 1];
P1 = [1 0 0 0; 0 1 0 0; 0 0 1 0];

det(U)
P2a = [U*W*V' (U(:,3))]
P2b = [U*W*V' -(U(:,3))]
P2c = [U*W'*V' (U(:,3))]
P2d = [U*W'*V' -(U(:,3))]

nx = x;
nx{1} = K\x{1};
nx{2} = K\x{2};

[Ma, X3Da]=triangulate(nx, P1, P2a);
[Mb, X3Db]=triangulate(nx, P1, P2b);
[Mc, X3Dc]=triangulate(nx, P1, P2c);
[Md, X3Dd]=triangulate(nx, P1, P2d);

%%
figure('Name', 'P2a')
X3Da = pflat(X3Da);

plot3(X3Da(1,:), X3Da(2,:), X3Da(3,:), '.')
%axis equal
hold on
plotcams({P1, P2a})
%%
figure('Name', 'P2b')
X3Db = pflat(X3Db);


plot3(X3Db(1,:), X3Db(2,:), X3Db(3,:), '.')
%axis equal
hold on
plotcams({P1, P2b})
%%
figure('Name', 'P2c')
X3Dc = pflat(X3Dc);

plot3(X3Dc(1,:), X3Dc(2,:), X3Dc(3,:), '.')
%axis equal
hold on
plotcams({P1, P2c})
%%
figure('Name', 'P2d')
X3Dd = pflat(X3Dd);

plot3(X3Dd(1,:), X3Dd(2,:), X3Dd(3,:), '.')
%axis equal
hold on
plotcams({P1, P2d})
%%
P2an = K*P2a;
xa = P2an * X3Da;
xa = pflat(xa);
im = imread('kronan2.JPG');
imagesc(im);
hold on
plot(xa(1,:),xa(2,:),'bo');
hold on
plot(x{2}(1,:),x{2}(2,:),'r*');
hold on
%%
P2bn = K*P2b;
xa = P2bn * X3Db;
xa = pflat(xa);
im = imread('kronan2.JPG');
imagesc(im);
hold on
plot(xa(1,:),xa(2,:),'bo');
hold on
plot(x{2}(1,:),x{2}(2,:),'r*');
hold on
%%
P2cn = K*P2c;
xa = P2cn * X3Dc;
xa = pflat(xa);
im = imread('kronan2.JPG');
imagesc(im);
hold on
plot(xa(1,:),xa(2,:),'bo');
hold on
plot(x{2}(1,:),x{2}(2,:),'r*');
hold on
%%
P2dn = K*P2d;
xa = P2dn * X3Dd;
xa = pflat(xa);
im = imread('kronan2.JPG');
imagesc(im);
hold on
plot(xa(1,:),xa(2,:),'bo');
hold on
plot(x{2}(1,:),x{2}(2,:),'r*');
hold on
%%
P1n = K*P1;
save('compEx1dataFrom', 'P1n','P2bn','x','X3Db');
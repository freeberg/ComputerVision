clear
load FundalMatrix.mat
load compEx1data.mat

P1 = [1 0 0 0; 0 1 0 0; 0 0 1 0];
e2 = null(F')
e2x = [0 -e2(3) e2(2); e2(3) 0 -e2(1); -e2(2) e2(1) 0]
P2 = [e2x*F e2]
%%
nx = x;
nx{1} = N1*x{1};
nx{2} = N2*x{2};
P1n = N1*P1
P2n = N1*P2

[M, X3D] = triangulate(nx, P1n, P2n);

xproj1 = pflat(P1n*X3D);
xproj2 = pflat(P2n*X3D);

nxproj1 = inv(N1)*xproj1;
nxproj2 = inv(N2)*xproj2;
%%
im = imread('kronan1.JPG');
imagesc(im)
hold on

plot(nxproj1(1,:), nxproj1(2,:), 'r*')
axis equal
plot(x{1}(1,:), x{1}(2,:), 'o')
%%
figure
im = imread('kronan2.JPG');
imagesc(im)
hold on

plot(nxproj2(1,:), nxproj2(2,:), 'r*')
axis equal
hold on
plot(x{2}(1,:), x{2}(2,:), 'o')
%%
figure
X3D = pflat(X3D)
plot3(X3D(1,:), X3D(2,:), X3D(3,:), '*')
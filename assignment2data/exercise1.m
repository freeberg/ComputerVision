clear
load compEx1data.mat
%%
figure
axis equal
plot3(X(1,:),X(2,:),X(3,:),'.')
axis equal
hold on
axis equal
plotcams(P)
axis equal
%%
%plot for camera 1
figure 
im = imread(imfiles{1});
imagesc(im);
hold on
visible = isfinite(x{1}(1,:));
axis equal
plot(x{3}(1,visible),x{1}(2,visible),'*','Markersize',2);
projpoints = P{1} * X
projpoints = pflat(projpoints);
axis equal
plot(projpoints(1,visible), projpoints(2,visible),'ro','Markersize',2)
%%
clear
load compEx1data.mat
%Transform matrix
T1 = [1 0 0 0; 0 4 0 0; 0 0 1 0; 1/10 1/10 0 1];
T2 = [1 0 0 0; 0 1 0 0; 0 0 1 0; 1/16 1/16 0 1];
%Inverse of Ts for the Ps
T1inv = inv(T1);
T2inv = (inv(T2));

%plot for T1
Xi1 = T1 * X;
Xi1 = pflat(Xi1);
axis equal
plot3(Xi1(1,:),Xi1(2,:),Xi1(3,:),'.','Markersize',2);
axis equal
P1=P;
P2=P;
for i=1:size(P,2)
    P1{i} = P{i} * T1inv;
end
hold on
axis equal
plotcams(P1)
axis equal

%plot for T2
figure
Xi2 = T2 * X;
Xi2 = pflat(Xi2);
axis equal
plot3(Xi2(1,:),Xi2(2,:),Xi2(3,:),'.','Markersize',2);
axis equal
for i=1:size(P,2)
    P2{i} = P{i} * T2inv;
end
hold on
axis equal
plotcams(P2)
axis equal


%plot for one camera

figure
im = imread(imfiles{3});
imagesc(im);
hold on
visible = isfinite(x{3}(1,:));
axis equal
plot(x{3}(1,visible),x{3}(2,visible),'*','Markersize',2);
projpoints = P2{3} * Xi2;
projpoints = pflat(projpoints);
plot(projpoints(1,visible), projpoints(2,visible),'ro','Markersize',2)

[r1 q1] = rq(P1{5})
[r2 q2] = rq(P2{5})

K1 = r1;
K2 = r2;

K1=K1./K1(3,3)
K2=K2./K2(3,3)

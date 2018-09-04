clear
load compEx3data.mat

%Camera 1
xy1mean = mean(x{1}(1:2,:),2);
xy1std = std(x{1}(1:2,:),0,2);

N = [1/xy1std(1) 0 -1/xy1std(1)*xy1mean(1); 0 1/xy1std(2) -1/xy1std(2)*xy1mean(2); 0 0 1];
nx = N*x{1};
pflat(nx)
axis equal
plot(nx(1,:),nx(2,:),'.')

[n, kx] = size(Xmodel)
[m, ky] = size(x{1})

M=[];
for col=1:kx
    for row=1:n
        M=[M; zeros(1,4*(row-1)),transpose(Xmodel(:,col)),1,zeros(1,4*(3-row)),zeros(1,col-1),-nx(row,col),zeros(1,(kx-col))];
    end
end

[U,S,V] = svd(M);
eigen = transpose(S)*S;
mineigen = min(eigen(eigen>0))

solu = V(:,end);
norm(M*solu)

P1 = [transpose(solu(1:4,1));transpose(solu(5:8,1));transpose(solu(9:12,1))];
%place projection infront of camera
P1 = -1*P1;
P1 = inv(N)*P1;

X = [Xmodel;ones(1,37)]
x1ImPoints = P1 * X;
x1ImPoints = pflat(x1ImPoints);
figure
im = imread('cube1.JPG');
imagesc(im);
hold on
plot(x{1}(1,:),x{1}(2,:),'o')
plot(x1ImPoints(1,:), x1ImPoints(2,:),'*');


%%
load compEx3data.mat

%Camera 2
xy1mean = mean(x{2}(1:2,:),2);
xy1std = std(x{2}(1:2,:),0,2);

N = [1/xy1std(1) 0 -1/xy1std(1)*xy1mean(1); 0 1/xy1std(2) -1/xy1std(2)*xy1mean(2); 0 0 1];
nx = N*x{2};
pflat(nx)
axis equal

[n, kx] = size(Xmodel)
[m, ky] = size(x{1})

M=[];
for col=1:kx
    for row=1:n
        M=[M; zeros(1,4*(row-1)),transpose(Xmodel(:,col)),1,zeros(1,4*(3-row)), ...
            zeros(1,col-1),-nx(row,col),zeros(1,(kx-col))];
    end
end

[U,S,V] = svd(M);
eigen = transpose(S)*S;
mineigen = min(eigen(eigen>0))

solu = V(:,end);
norm(M*solu)

P2 = [transpose(solu(1:4,1));transpose(solu(5:8,1));transpose(solu(9:12,1))];
%place projection infront of camera
P2 = -1*P2;
P2 = inv(N)*P2;

X = [Xmodel;ones(1,37)]
x1ImPoints = P2 * X;
x1ImPoints = pflat(x1ImPoints);
figure
im = imread('cube2.JPG');
imagesc(im);
hold on
plot(x{2}(1,:),x{2}(2,:),'o')
plot(x1ImPoints(1,:), x1ImPoints(2,:),'*');
%%
plot3([Xmodel(1,startind );  Xmodel(1,endind )],...
      [Xmodel(2,startind );  Xmodel(2,endind )],...
      [Xmodel(3,startind );  Xmodel(3,endind)],'b-');
hold on
plot3(Xmodel(1,:),Xmodel(2,:),Xmodel(3,:),'x')
cams = {1,2};
cams{1} = P1;
cams{2} = P2;
plotcams(cams)

%plot for camera centers
cc1 = null(P1);
cc2 = null(P2);
cc1 = pflat(cc1);
cc2 = pflat(cc2);
plot3(cc1(1,:),cc1(2,:),cc1(3,:),'g*')
plot3(cc2(1,:),cc2(2,:),cc2(3,:),'b*')

%%
%Inner parameters for Camera 1
[K Q] = rq(P1)

%%
%Save camera matrixes for exercise5
save('cameraP1.mat','P1');
save('cameraP2.mat','P2');
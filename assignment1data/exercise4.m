figure
image1 = imread('compEx4im1.JPG');
imagesc(image1);

figure
image2 = imread('compEx4im2.JPG');
imagesc(image2);

load compEx4.mat
camcent1 = null(P1)
camcent2 = null(P2)

%%

X = [0; 0; 1] 

u3D = pflat(U);
camcent1 = pflat(camcent1);
camcent2 = pflat(camcent2);
camcents = [camcent1, camcent2]
prinax1 = P1(3,1:3);
prinax2 = P2(3,1:3);
prinax1 = prinax1/norm(prinax1)
prinax2 = prinax2/norm(prinax2)
figure
plot3(u3D(1,:),u3D(2,:),u3D(3,:),'.','Markersize',2)
hold on
plot3(camcent1(1), camcent1(2), camcent1(3),'x')
plot3(camcent2(1), camcent2(2), camcent2(3),'x')
quiver3(camcent1(1), camcent1(2), camcent1(3), prinax1(1), prinax1(2), prinax1(3),1)
quiver3(camcent2(1), camcent2(2), camcent2(3), prinax2(1), prinax2(2), prinax2(3),1)

%%
figure
image1 = imread('compEx4im1.JPG');
imagesc(image1);
hold on
impoints = P1*U;
impoints = pflat(impoints)
plot3(impoints(1,:),impoints(2,:),impoints(3,:),'x')

figure
image2 = imread('compEx4im2.JPG');
imagesc(image2);
hold on
impoints2 = P2*U;
impoints2 = pflat(impoints2)
plot3(impoints2(1,:),impoints2(2,:),impoints2(3,:),'x')
%%

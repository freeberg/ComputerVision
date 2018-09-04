function disp_result(im0,P,segm0,depths,tol,sc)
[K,R] = rq(P);
%rescale
im0 = imresize(im0,sc);
segm0 = imresize(segm0,sc);
K(1:2,:)= sc*K(1:2,:);
P = K*R;

[xx,yy] = meshgrid(1:size(im0,2),1:size(im0,1));
impoints = inv(K)*pextend([xx(segm0)'; yy(segm0)']);
depths = depths(segm0)';

%Skapa triangulering.
nodes = zeros(size(im0(:,:,1)));
nodes(segm0) = 1:length(nodes(segm0));

corner1 = [nodes(1:end-1,1:end-1)]; 
corner2 = [nodes(1:end-1,2:end)]; 
corner3 = [nodes(2:end,1:end-1)]; 

keep_ind = (corner1 ~= 0) & (corner2 ~= 0) & (corner3 ~= 0);

corner1 = corner1(keep_ind);
corner2 = corner2(keep_ind);
corner3 = corner3(keep_ind);

depth1 = depths(corner1);
depth2 = depths(corner2);

cutedges = abs(depth1-depth2) > tol;

depth1 = depths(corner3);
depth2 = depths(corner2);

cutedges = cutedges | abs(depth1-depth2) > tol;

depth1 = depths(corner3);
depth2 = depths(corner1);

cutedges = cutedges | abs(depth1-depth2) > tol;

tri = [corner1(~cutedges) corner2(~cutedges) corner3(~cutedges)];
imgray = rgb2gray(im0);

color = imgray(1:end-1,1:end-1);
color = color(keep_ind);
tricolor = color(~cutedges);


corner1 = [nodes(1:end-1,2:end)];
corner3 = [nodes(2:end,1:end-1)];
corner2 = [nodes(2:end,2:end)];

keep_ind = (corner1 ~= 0) & (corner2 ~= 0) & (corner3 ~= 0);

corner1 = corner1(keep_ind);
corner2 = corner2(keep_ind);
corner3 = corner3(keep_ind);


depth1 = depths(corner1);
depth2 = depths(corner2);

cutedges = abs(depth1-depth2) > tol;

depth1 = depths(corner3);
depth2 = depths(corner2);

cutedges = cutedges | abs(depth1-depth2) > tol;

depth1 = depths(corner3);
depth2 = depths(corner1);

cutedges = cutedges | abs(depth1-depth2) > tol;

color = imgray(1:end-1,1:end-1);
color = color(keep_ind);
tricolor = [tricolor; color(~cutedges)];

tri = [tri; [corner1(~cutedges) corner2(~cutedges) corner3(~cutedges)]];

points3D = impoints.*repmat(depths,[3 1]);

figure;
trisurf(tri,points3D(1,:),points3D(2,:),points3D(3,:),tricolor,'EdgeColor','None');
colormap gray;
hold on;

%draw camera
lambda = 1;
[xx,yy] = meshgrid([1:(size(im0,2)+1)]-0.5,[1:(size(im0,1)+1)]-0.5);
corners = pextend([-0.5 -0.5 size(im0,2)+0.5 size(im0,2)+0.5;...
    -0.5 size(im0,1)+0.5 size(im0,1)+0.5 -0.5]);
imgray = rgb2gray(im0);
[K,R] = rq(P);
points1 = R(:,1:3)'*(lambda*inv(K)*pextend([xx(:)';yy(:)'])-repmat(R(:,4),[1 length(xx(:))]));
corners1 = R(:,1:3)'*(lambda*inv(K)*corners-repmat(R(:,4),[1 4]));
nodenr = zeros(size(xx));
nodenr(:) = 1:length(nodenr(:));
tri1 = [];
color1 = [];
corn1 = nodenr(1:end-1,1:end-1);
corn2 = nodenr(1:end-1,2:end);
corn3 = nodenr(2:end,1:end-1);
tri1 = [tri1; corn1(:) corn2(:) corn3(:)];
color1 = [color1; imgray(:)];
corn1 = nodenr(2:end,2:end);
corn2 = nodenr(1:end-1,2:end);
corn3 = nodenr(2:end,1:end-1);
tri1 = [tri1; corn1(:) corn2(:) corn3(:)];
color1 = [color1; imgray(:)];

trisurf(tri1,points1(1,:),points1(2,:),points1(3,:),color1,'edgecolor','none');

c1 = -R(:,1:3)'*R(:,4);
plot3(c1(1),c1(2),c1(3),'*');
plot3(corners1(1,[1:4 1]),corners1(2,[1:4 1]),corners1(3,[1:4 1]),'-','linewidth',2);
plot3([corners1(1,:)' ones(4,1)*c1(1)]',[corners1(2,:)' ones(4,1)*c1(2)]',[corners1(3,:)' ones(4,1)*c1(3)]','b-');
axis equal
hold off;

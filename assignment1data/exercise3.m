load compEx3.mat


figure
plot([startpoints(1,:);endpoints(1,:)] ...
    ,[startpoints(2,:);endpoints(2,:)], '-b')


H1 = [sqrt(3), -1, 1; 1, sqrt(3), 1; 0, 0, 2]
H2 = [1, -1, 1; 1, 1, 0; 0, 0, 1]
H3 = [1, 1, 0; 0, 2, 0; 0, 0, 1]
H4 = [sqrt(3), -1, 1; 1, sqrt(3), 1; 1/4, 1/2, 2]

startpoints = [startpoints; ones(1,size(startpoints,2))]
endpoints = [endpoints; ones(1,size(endpoints,2))]

figure
plotgrid(H1*startpoints, H1*endpoints)
figure
plotgrid(H2*startpoints, H2*endpoints)
figure
plotgrid(H3*startpoints, H3*endpoints)
figure
plotgrid(H4*startpoints, H4*endpoints)

function plotgrid(sp, ep)
axis equal
sp=pflat(sp);
ep=pflat(ep);
    plot([sp(1,:);ep(1,:)] ...
        , [sp(2,:);ep(2,:)], '-b')
end
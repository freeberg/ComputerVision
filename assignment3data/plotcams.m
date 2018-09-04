function plotcams(P)

c = zeros(4,length(P));
v = zeros(3,length(P));
for i = 1:length(P);
    c(:,i) = null(P{i});
    v(:,i) = P{i}(3,1:3);
end
c = c./repmat(c(4,:),[4 1]);
quiver3(c(1,:),c(2,:),c(3,:),v(1,:), v(2,:), v(3,:),['r','-'],'LineWidth',1.5,'MaxHeadSize',1.5);
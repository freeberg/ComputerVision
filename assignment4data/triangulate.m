function [M, X3D]=triangulate(xt, P1, P2)
clear X3D;
M={};

[n, kx] = size(xt{1});
eqnbr = 0;
for i=1:kx
    eqnbr = eqnbr+1;
    M{i}=[P1, -xt{1}(:,i), zeros(3, 1); P2, zeros(3, 1), -xt{2}(:,i)];
    
end
X3D = [];
for j=1:kx
    [UT,ST,VT] = svd(M{j});
    solu = VT(1:4,end);
    X3D = [X3D, solu(1:4,1)];
end

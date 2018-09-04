clear
load compEx1data.mat

X = pflat(X);
meanX = mean(X, 2);

Xtilde = (X - repmat(meanX ,[1 size(X ,2)]));
M = Xtilde (1:3 ,:)* Xtilde (1:3,:)';

[V,D] = eig(M) %V eigenvectors, D eigenvalues in diag
[minEig, I] = min(diag(D));
minEigVec = V(:,I)

d = -(minEigVec(1)*meanX(1) + minEigVec(2)*meanX(2) + minEigVec(3)*meanX(3))

tlsPlane = [minEigVec; d]
tlsPlane = tlsPlane./norm(tlsPlane (1:3))

RMStls = sqrt(sum((tlsPlane' * X).^2)/ size(X,2))
hist(abs(tlsPlane' * X), 100);
%%
[m, n] = size(X);
bestInl = 0;
bestplane = []
bestInliers = [];
for i=1:10 
   perm = randperm(n, 3);
   plane = null(X(:, perm)');
   plane = plane./norm(plane(1:3));
   inliers = abs(plane'*X) <= 0.1;
   amtInl = sum(inliers(:) == 1);
   if (amtInl > bestInl)
       bestplane = plane;
       bestInl = amtInl;
       bestInliers = inliers;
   end
end
figure
RMSransac = sqrt(sum((bestplane' * X).^2)/ size(X,2))
hist(abs(bestplane' * X), 100);
%%
Xinl = X(:, bestInliers);

meanXinl = mean(Xinl, 2);

Xinltilde = (Xinl - repmat(meanXinl ,[1 size(Xinl ,2)]));
M = Xinltilde(1:3 ,:) * Xinltilde(1:3,:)';

[V,D] = eig(M) %V eigenvectors, D eigenvalues in diag
[minEig, I] = min(diag(D));
minEigVec = V(:,I)

d = -(minEigVec(1)*meanXinl(1) + minEigVec(2)*meanXinl(2) + minEigVec(3)*meanXinl(3))
pXinl = pflat(Xinl);
tlsPlane = [minEigVec; d]
tlsPlane = tlsPlane./norm(tlsPlane (1:3))

RMSinltls = sqrt(sum((tlsPlane' * pXinl).^2)/ size(Xinl,2))
hist(abs(tlsPlane' * pXinl), 100);
%%
figure
im = imread('house1.jpg');
imagesc(im);
hold on
xproj = P{1} * pXinl;
xproj = pflat(xproj);
plot(xproj(1,:), xproj(2,:), '*')

figure
im = imread('house2.jpg');
imagesc(im);
hold on
xproj = P{2} * pXinl;
xproj = pflat(xproj);
plot(xproj(1,:), xproj(2,:), '*')

%%
P2n = inv(K) * P{2}
xk = inv(K) * x;

tlsPlane = pflat(tlsPlane);
H = (P2n(:,1:3) - P2n(:,4)*tlsPlane(1:3)');

xh = H * xk
xh = pflat(K * xh)


figure
im = imread('house1.jpg');
imagesc(im);
hold on;
plot(x(1,:), x(2,:), 'r*')

figure
im = imread('house2.jpg');
imagesc(im);
hold on;
plot(xh(1,:), xh(2,:), 'r*')
clear
load compEx3data.mat
load compEx1data.mat

nx1 = K\x{1};
nx2 = K\x{2};

[m, kx] = size(nx1);

M=[];
%Create the 8-point M-matrix
for col=1:kx
    xx = nx2(:,col)*nx1(:,col)';
    M(col,:) = xx(:)';
end
M(1,:)
[U,S,V] = svd(M);

eigen = transpose(S)*S;
mineigen = min(eigen(eigen>0))
solu = V(:,end);
norm = norm(M*solu)

E = reshape(solu, [3 3]);
detE = det(E)
[U,S,V] = svd(E);
if det(U*V')>0
E = U*diag ([1 1 0])*V';
else
V = -V;
E = U*diag ([1 1 0])*V';
end
save('EssenMatrix.mat','U','V');

detE = det(E)
plot(diag(nx2'*E*nx1))
%E = E./E(3,3)
%%

F = K'\E/K; % (K^T)⁻1 * E * K^-1
plot(diag(x{2}'*F*x{1}));
l = F*x{1};
l = l./sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2 ,[3 1]));
figure
perm = randperm(20);
im = imread('kronan2.JPG');
imagesc(im)
hold on
rital(l(:,perm(1:20)),'b-')
plot(x{2}(1,perm(1:20)), x{2}(2,perm(1:20)),'r*')

%%
hist(abs(sum(l.*x{2})), 100);
meanDistance = mean(abs(sum(l.*x{2})))


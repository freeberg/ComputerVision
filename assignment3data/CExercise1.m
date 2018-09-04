clear
load compEx1data.mat

%calculate the image points mean and std
xy1mean = mean(x{1}(1:2,:),2);
xy1std = std(x{1}(1:2,:),0,2);
xy2mean = mean(x{2}(1:2,:),2);
xy2std = std(x{2}(1:2,:),0,2);

%find N matrixes
N1 = [1/xy1std(1) 0 -1/xy1std(1)*xy1mean(1); 0 1/xy1std(2) -1/xy1std(2)*xy1mean(2); 0 0 1];
N2 = [1/xy2std(1) 0 -1/xy2std(1)*xy2mean(1); 0 1/xy2std(2) -1/xy2std(2)*xy2mean(2); 0 0 1];
%normalize the image points
nx1 = N1*x{1};
nx2 = N2*x{2};

[m, kx] = size(x{1});

M=[];
%Create the 8-point M-matrix
for col=1:kx
    xx = nx2(:,col)*nx1(:,col)';
    M(col,:) = xx(:)';
end
M(1,:)
[U,S,V] = svd(M);



%since not zero we ensure that it is
solu = V(:,end);
Fn = reshape(solu, [3 3])
detFnBefore = det(Fn)
[U,S,V] = svd(Fn);
S(3,end)=0;
Fn = U*S*transpose(V)

%Check so min singular value and ||Mv|| is close to zero
eigen = transpose(S)*S;
mineigen = min(eigen(eigen>0))
solu = V(:,end);
normFn = norm(Fn*solu)

%Check so x2T*F*x1 is roughly zero
detFnAfter = det(Fn)
plot(diag(nx2'*Fn*nx1));

%Compute unnormalized F-matrix using formula from exer3
F = transpose(N2)*Fn*N1;
plot(diag(x{2}'*F*x{1}));
l = F*x{1};
l = l./sqrt(repmat(l(1 ,:).^2 + l(2 ,:).^2 ,[3 1]));
%%
perm = randperm(20);
im = imread('kronan2.JPG');
imagesc(im)
hold on
rital(l(:,perm(1:20)),'b-')
plot(x{2}(1,perm(1:20)), x{2}(2,perm(1:20)),'r*')
%%
hist(abs(sum(l.*x{2})) ,100);
meanDistance = mean(abs(sum(l.*x{2})))
F = F./F(3,3)
save('FundalMatrix.mat','F','N1','N2');
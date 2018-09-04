clear
load compEx1dataFrom.mat

u = x;
U = X3Db;
im1 = imread('kronan1.JPG');
im2 = imread('kronan2.JPG');

P = {P1n, P2bn};

[err ,res] = ComputeReprojectionError(P,U,u);
figure
hist(res, 100)

lambda = 10000;

%%

n = 10;
ii = 1:n;


%%

for j=1:12
    P = {P1n, P2bn};
    u = x;
    U = X3Db;
    errors = zeros(1, n + 1); 
    [errors(1), res] = ComputeReprojectionError(P,U,u);
    for i=ii
    [r,J] = LinearizeReprojErr(P,U,u);
    C = J'*J+lambda*speye(size(J,2));
    c = J'*r;
    deltav = -C\c;
    
    [P, U] = update_solution(deltav,P,U);
    
    [errors(i+1), res] = ComputeReprojectionError(P,U,u);
    end
    lambda = lambda / 100;
    plot(0:n,errors);
    hold on
end

%%
figure
hist(res, 100)
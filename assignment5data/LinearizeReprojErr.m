function [r,J] = LinearizeReprojErr(P,U,u)
%Linearizes the reprojection error functions for the calibrated structure
%from motion problem.
%
%Inputs:    P - cell containing the current estimates of the cameras.
%           U - 4xn matrix containing the current estimate of the
%           structure (3D-points)
%           u - cell containing the imagepoints visible in each image.
%
%Outputs:   r - values of all residuals at the current solution.
%           J - Jacobian at the current solution.


%Normalize cameras
for i=1:length(P);
    [K,R] = rq(P{i});
    u{i} = inv(K)*u{i};
    P{i} = R;
end

numpts = size(U,2);
%Bas f�r tangentplanet till rotationsm�ngfalden.
Ba = [0 1 0; -1 0 0; 0 0 0];
Bb = [0 0 1; 0 0 0; -1 0 0];
Bc = [0 0 0; 0 0 1; 0 -1 0];

da1 = cell(size(P));
db1 = cell(size(P));
dc1 = cell(size(P));
dt11 = cell(size(P));
dt21 = cell(size(P));
dt31 = cell(size(P));
dU11 = cell(size(P));
dU21 = cell(size(P));
dU31 = cell(size(P));
da2 = cell(size(P));
db2 = cell(size(P));
dc2 = cell(size(P));
dt12 = cell(size(P));
dt22 = cell(size(P));
dt32 = cell(size(P));
dU12 = cell(size(P));
dU22 = cell(size(P));
dU32 = cell(size(P));


U = pflat(U);
U = U(1:3,:);

for i=1:length(P);
    %ber�kna derivator f�r b�da residualeran i alla bilder
    %a,b,c - rotations parametrar f�r kamera i
    %U1,U2,U3 - 3d punkt parametrar
    %t1,t2,t3 - translations parametrar kameran.
    vis = isfinite(u{i}(1,:));
    R0 = P{i}(:,1:3);
    %   t0 = repmat(P{i}(:,4),[1 size(u{i},2)]);
    t0 = P{i}(:,4);
    da1{i} =  (Ba(1,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(1,:)*U(:,vis)+t0(1))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Ba(3,:)*R0*U(:,vis));

    da2{i} = (Ba(2,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(2,:)*U(:,vis)+t0(2))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Ba(3,:)*R0*U(:,vis));
    
    db1{i} = (Bb(1,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(1,:)*U(:,vis)+t0(1))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Bb(3,:)*R0*U(:,vis));
    
    db2{i} = (Bb(2,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(2,:)*U(:,vis)+t0(2))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Bb(3,:)*R0*U(:,vis));
    
    dc1{i} = (Bc(1,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(1,:)*U(:,vis)+t0(1))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Bc(3,:)*R0*U(:,vis));

    dc2{i} = (Bc(2,:)*R0*U(:,vis))./(R0(3,:)*U(:,vis)+t0(3)) - ...
        (R0(2,:)*U(:,vis)+t0(2))./((R0(3,:)*U(:,vis)+t0(3)).^2).*(Bc(3,:)*R0*U(:,vis));
    
    dU11{i} = R0(1,1)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(1,:)*U(:,vis) + t0(1))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,1);

    dU12{i} = R0(2,1)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(2,:)*U(:,vis) + t0(2))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,1);
    
    dU21{i} = R0(1,2)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(1,:)*U(:,vis) + t0(1))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,2);

    dU22{i} = R0(2,2)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(2,:)*U(:,vis) + t0(2))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,2);
    
    dU31{i} = R0(1,3)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(1,:)*U(:,vis) + t0(1))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,3);
    
    dU32{i} = R0(2,3)./(R0(3,:)*U(:,vis) + t0(3)) - ...
        (R0(2,:)*U(:,vis) + t0(2))./((R0(3,:)*U(:,vis) + t0(3)).^2).*R0(3,3);

    dt11{i} = 1./(R0(3,:)*U(:,vis)+t0(3));
    dt12{i} = zeros(size(dt11{i}));
    
    dt21{i} = zeros(size(dt11{i}));
    dt22{i} = 1./(R0(3,:)*U(:,vis)+t0(3));
    
    dt31{i} = -(R0(1,:)*U(:,vis)+t0(1))./((R0(3,:)*U(:,vis)+t0(3)).^2);
    dt32{i} = -(R0(2,:)*U(:,vis)+t0(2))./((R0(3,:)*U(:,vis)+t0(3)).^2);
        
end
row = [];
col = [];
data = [];
resnum = 0;
B = [];
%fprintf('\tSetting up system, camera: ');
%fprintf('\tSetting up system. ');
%Allokera hela vektorn f�rst tv� residualer f�r varje bild punkt
resnr = 0;
for i = 1:length(P);
    resnr = resnr + sum(isfinite(u{i}(1,:)));
end
row = zeros(2*resnr*(6+3),1);
col = zeros(2*resnr*(6+3),1);
data = zeros(2*resnr*(6+3),1);
B = zeros(2*resnr,1);
lastentry = 0;
lastentryB = 0;
for i = 1:length(P);
    %fprintf('%d ',i);
    uu = u{i};
    vis = find(isfinite(uu(1,:)));
    
    %F�rsta residualen:
    %3D-punkt parametrar:
    %U1-koeff
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [(vis-1)*3+1]';
    data(lastentry+[1:length(vis)]) = dU11{i}';
    lastentry = lastentry+length(vis);
    
    %U2-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; [(vis-1)*3+2]'];
    %data = [data; dU21{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [(vis-1)*3+2]';
    data(lastentry+[1:length(vis)]) = dU21{i}';
    lastentry = lastentry+length(vis);
    
    %U3-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; [vis*3]'];
    %data = [data; dU31{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [vis*3]';
    data(lastentry+[1:length(vis)]) = dU31{i}';
    lastentry = lastentry+length(vis);
    
    %Kameraparametrar
    %a-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+1)*ones(length(vis),1)];
    %data = [data; da1{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+1)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = da1{i}';
    lastentry = lastentry+length(vis);
    
    %b-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+2)*ones(length(vis),1)];
    %data = [data; db1{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+2)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = db1{i}';
    lastentry = lastentry+length(vis);
    
    %c-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+3)*ones(length(vis),1)];
    %data = [data; dc1{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+3)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dc1{i}';
    lastentry = lastentry+length(vis);
    
    %t_1-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+4)*ones(length(vis),1)];
    %data = [data; dt11{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+4)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt11{i}';
    lastentry = lastentry+length(vis);
    
    %t_2-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+5)*ones(length(vis),1)];
    %data = [data; dt21{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+5)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt21{i}';
    lastentry = lastentry+length(vis);
    
    %t_3-koeff
    %row = [row; resnum+[1:2:2*length(vis)]'];
    %col = [col; (3*numpts+i*6)*ones(length(vis),1)];
    %data = [data; dt31{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[1:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+i*6)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt31{i}';
    lastentry = lastentry+length(vis);
    
    %Andra residualen:
    %3D-punkt parametrar:
    %U1-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; [(vis-1)*3+1]'];
    %data = [data; dU12{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [(vis-1)*3+1]';
    data(lastentry+[1:length(vis)]) = dU12{i}';
    lastentry = lastentry+length(vis);
    
    
    %U2-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; [(vis-1)*3+2]'];
    %data = [data; dU22{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [(vis-1)*3+2]';
    data(lastentry+[1:length(vis)]) = dU22{i}';
    lastentry = lastentry+length(vis);
    
    %U3-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; [vis*3]'];
    %data = [data; dU32{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = [vis*3]';
    data(lastentry+[1:length(vis)]) = dU32{i}';
    lastentry = lastentry+length(vis);

    %Kameraparametrar
    %a-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+1)*ones(length(vis),1)];
    %data = [data; da2{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+1)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = da2{i}';
    lastentry = lastentry+length(vis);
    
    %b-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+2)*ones(length(vis),1)];
    %data = [data; db2{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+2)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = db2{i}';
    lastentry = lastentry+length(vis);
    
    %c-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+3)*ones(length(vis),1)];
    %data = [data; dc2{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+3)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dc2{i}';
    lastentry = lastentry+length(vis);
    
    %t_1-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+4)*ones(length(vis),1)];
    %data = [data; dt12{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+4)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt12{i}';
    lastentry = lastentry+length(vis);
    
    %t_2-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+(i-1)*6+5)*ones(length(vis),1)];
    %data = [data; dt22{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+(i-1)*6+5)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt22{i}';
    lastentry = lastentry+length(vis);
        
    %t_3-koeff
    %row = [row; resnum+[2:2:2*length(vis)]'];
    %col = [col; (3*numpts+i*6)*ones(length(vis),1)];
    %data = [data; dt32{i}'];
    row(lastentry+[1:length(vis)]) = resnum+[2:2:2*length(vis)]';
    col(lastentry+[1:length(vis)]) = (3*numpts+i*6)*ones(length(vis),1);
    data(lastentry+[1:length(vis)]) = dt32{i}';
    lastentry = lastentry+length(vis);
        
    resnum = resnum+2*length(vis);
    
    %Konstant-termer
    btmp = zeros(2*length(vis),1);
    %F�rsta residualen
    btmp(1:2:end) = (P{i}(1,:)*pextend(U(:,vis)))./(P{i}(3,:)*pextend(U(:,vis)))-uu(1,vis);
    %Andra residualen
    btmp(2:2:end) = (P{i}(2,:)*pextend(U(:,vis)))./(P{i}(3,:)*pextend(U(:,vis)))-uu(2,vis);
    %B = [B; btmp];
    B(lastentryB+[1:length(btmp)]) = btmp;
    lastentryB = lastentryB+length(btmp);
    
end

A = sparse(row,col,data);
%l�s koordinatsystemet
%F�rsta kameran konstant
%och f�rsta koordinaten i f�rsta punkten konstant
A = A(:,[1:3*numpts 3*numpts+7:end]);
A = A(:,[2:end]);

r = B;
J = A;
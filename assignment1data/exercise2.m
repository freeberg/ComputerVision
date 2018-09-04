image = imread('compEx2.JPG');
imagesc(image);
colormap gray;
hold on;
load compEx2.mat;

plotpoints(p1);
plotpoints(p2);
plotpoints(p3);

linje1=null(transpose(p1))
linje2=null(transpose(p2))
linje3=null(transpose(p3))


rital([linje1,linje2,linje3])

intsec=null([transpose(linje2);transpose(linje3)]);
intsec=pflat(intsec)
plotpoints(intsec);

t = abs(linje1(1)*intsec(1)+linje1(2)*intsec(2)+linje1(3))
n = sqrt(linje1(1)^2 + linje1(2)^2)
d = t/n


function plotpoints(p)
   plot(p(1,:),p(2,:),'x');
   set(findall(gca, 'Type', 'Line'),'LineWidth',1.5); 
end
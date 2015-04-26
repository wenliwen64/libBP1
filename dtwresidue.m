function [residue velocity warpsyn]=dtwresidue(syn,data,eps,del,flag)

d1=(syn-mean(syn))/std(syn);
d2=(data-mean(data))/std(data);
n=length(syn);
discut=1:n;

if nargin == 2 
eps=0.5;
del=0.1;
flag=0;
end

if nargin == 4 

flag=0;
end



epsilon =eps*min(std(d1), std(d2));
delta =round( del*n);


   [a b sim]=lcsMatching(d2', d1', delta,epsilon, 4);
      
    
   
d1in=interp1(discut(a),d1(b),discut,'linear','extrap');
discut1=interp1(discut(b),discut(a),discut,'linear','extrap');
discut1=smooth(discut1,50,'moving');
% k=round(n/(max(discut)-min(discut))*20);
k=round(n/(max(discut)-min(discut))*20);

for i=1:n-k
    ddis(i)=(discut1(i+k)-discut1(i))/(discut(i+k)-discut(i));
end

velocity=ddis;
% residue=d2-d1in;
residue=resample(d2(a)-d1(b),length(d2),length(a));
% residue(end-round(end/20):end)=0;
% residue(1:round(end/20))=0;

warpsyn=d1in;
% figure(200);
% plot(discut(1:n-k),ddis);
% plot(discut(1:n-k),smooth(ddis,5,'sgolay'));
if flag==1
figure(21);
set(gcf,'Position',[500 600 1600 600]);
subplot(411);
plot(discut,d1,discut,d2,'r');
subplot(412);
plot(discut,d1in,discut,d2,'r',discut,d2-d1in,'g');
subplot(413);
plot(discut(1:n-k),ddis);
subplot(414);
plot(discut,residue);
title(['residue' num2str(std(residue)) 'data' num2str(std(d2)) 'warpsyn' num2str(std(d1in))]);
figure(23);
plot(a,'r.');
hold on;
plot(b,'b.');
figure(25);
subplot(311);
plot(discut(a),d2(a),'r',discut(b),d1(b),'b',discut(a),d1(b),'c',discut(b),d2(a),'c');
subplot(312);
plot(discut,d2,'r',discut,d1,'b');
subplot(313);
plot(1:length(a),d2(a),'r',1:length(b),d1(b),'b');
% figure(23);

end
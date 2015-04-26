function plotAlllocal(DataFile,scale)
close all;
load(DataFile);
x=ret.xori;
sr=ret.sr;
fl=ret.fl;
fh=ret.fh;
sarr=ret.sarr;
[n m]=size(x);
[BB,AA]=butter(4,[fl fh]/(sr/2));
for i=1:n  
    x(i,:)=filter(BB,AA,x(i,:));
end
for i=1:n 
        x(i,:)=x(i,:)/std(x(i,sarr*ret.sr:(sarr+3)*ret.sr));
end
[n m]=size(x);
t=(1:m)/sr;
figure(30);
hold on;
for i=1:n
          plot(t,x(i,1:m)*scale+i-1,'b');
end

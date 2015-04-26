b0=100;
e0=170;
figure(7);
ttt=(1:length(ret.x))/ret.sr;
% plot(ttt,ret1.x(1,:),ttt,ret2.x(1,:),'r');
tmp=zeros(1,1);
len=length(ret.x(1,b0*ret.sr:e0*ret.sr));
for i=1:nel
    tmp(i,1:len)=ret2.x(i,b0*ret.sr:e0*ret.sr);
    tmpf=fft(tmp(i,:));
    tmpa=abs(tmpf);
    tmpph=tmpf./tmpa;
  
    r=2*pi*rand(1,(length(tmpph)-1)/2);
    tmpph(2:(length(tmpph)-1)/2+1)=cos(r)+1i*sin(r);
    tmpph((length(tmpph)-1)/2+2:length(tmpph))=fliplr(tmpph(2:(length(tmpph)-1)/2+1))';
    tmp(i,1:len)=ifft(tmpph.*tmpa);
    
end
ttt1=ttt(b0*ret.sr:e0*ret.sr);
figure(7);
plot(ttt,ret1.x(1,:),ttt1,tmp(1,:),'r');
xlim([b0 e0]);
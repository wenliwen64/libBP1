clear all;
sr = 200; %sampling rate /sec
%t = linspace(0,1,sr);
%x = sin(2*pi*45*t)+sin(2*pi*55*t)+sin(2*pi*60*t);% construct the signal
dur=10;
nel=12;
x1=zeros(nel,dur*sr);
k=1;
t=linspace(0,dur,dur*sr);
filename='2821044o4.p';
[xtext pr]=load_bbdata([filename '13' ],[],dur);
start=pr.t0;
for i=1:13
    
    if i==4
        continue;
    end
    if i<=9
   %[x(k,:) pr]=load_bbdata(['2721715j1.p0' num2str(i)],[],dur);
    %[x(k,:) pr]=load_bbdata(['3192141g4.p0' num2str(i)],[],dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p0' num2str(i)],[2004 276 07 04 28 575 0],dur);
    [x1(k,:) pr]=load_bbdata([filename '0' num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p0' num2str(i)],[],dur);
    else
  % [x(k,:) pr]=load_bbdata(['2721715j1.p' num2str(i)],[],dur);  
    %[x(k,:) pr]=load_bbdata(['2710200e2.p' num2str(i)],[],dur);
    [x1(k,:) pr]=load_bbdata([filename num2str(i)],start,dur);
    %[x(k,:) pr]=load_bbdata(['2760704k4.p' num2str(i)],[2004 276 07 04 28 575 0],dur);
    %[x(k,:) pr]=load_bbdata(['2910557b4.p' num2str(i)],[],dur);
    end
    x1(k,:)=x1(k,:)-mean(x1(k,:));
    k=k+1;
%     if i==13
%         start=pr.t0;
%     end
end

wd=2;
tl=2.8;
th=tl+wd;
x=zeros(nel,wd*sr);
for i=1:nel
x(i,1:wd*sr)=x1(i,wd*sr:4*sr-1);
end

for i=1:nel
    X(i,:)=fft(x(i,:));
    
end

 s1=linspace(0,sr,wd*sr);
          [s c ph ci phi]=cmtm1(x(1,:),x(2,:),0.005,1.5,0,0,0);
      figure(14);
      subplot(3,2,1);
      plot(s,abs(c),'b-');
       xlim([0 40]);
       subplot(3,2,2);
      plot(s,ph,'r-');
           
      subplot(3,2,3);
      plot(s1,abs(X(1,:).*X(2,:)),'b-');
      xlim([0 40]);
       subplot(3,2,4);
      plot(s1,angle(X(1,:).*X(2,:)),'r-');
      xlim([0 100]);
       subplot(3,2,5);
       [Coh,Pha,DPha] = cmtm2(x(1,:),x(2,:),1.5);
       plot(s1,abs(Coh),'b-');
       xlim([0 40]);
       subplot(3,2,6);
      plot(s1,Pha,'b-');
      xlim([0 100]);
      
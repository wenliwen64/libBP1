function [p1 p2 p1l p1r p2l p2r]=drawmax(Qm,percentdrop)
% plotting finding the percent drop contour for a colormap

% load /home/lsmeng/matlab/haitiUS/paper/distanceVsseparation_beam3to1/rup11.mat
% Qm=rupall;
% figure(20);
% hold on;
% pcolor(rupall);
% shading flat;
% percentdrop=0.1;
[n m]=size(Qm);
mindist=5;
 for i=1:n
     pi=ones(1,2);
     tmp=Qm(i,:);
     [value peaks]=findpeaks(tmp+0.001*max(tmp)*rand(1,length(tmp)),'minpeakdistance',mindist);
     np=length(peaks);
     tmp1=sortrows([peaks' tmp(peaks)'],2);
     peaks0=tmp1(:,1);
     low=min(tmp);
     pi(1)=peaks0(np);
     
     if np>=2%&&
         if abs(peaks0(np)-peaks0(np-1))>mindist
             M=2;
             
             pi(2)=peaks0(np-1);
             if pi(1)>pi(2)
                 tt=pi(1);
                 pi(1)=pi(2);
                 pi(2)=tt;
             end
             
         else
             M=1;
         end
     else
         M=1;
     end
 
     
    

     pil=ones(1,2);
     pir=ones(1,2);
     for k=1:M
         clear pr thes
         thes=(1-percentdrop)*(tmp(pi(k))-low)+low;
         count=1;
         for j=1:m-1
             
             if (tmp(j)-thes)*(tmp(j+1)-thes)<=0
                 pr(count)=j;
                 count=count+1;
             end
         end
         
         for j=1:length(pr)-1
             if (pr(j)-pi(k))*(pr(j+1)-pi(k))<=0
                 pil(k)=pr(j);
                 pir(k)=pr(j+1);
             end
         end
     end
     p1l(i)=pil(1);
     p1r(i)=pir(1);
     p2l(i)=pil(2);
     p2r(i)=pir(2);
     
     
%      peak2=localMaximum(-abs(tmp-thes));
 
 p1(i)=pi(1);
 p2(i)=pi(2);
 end
%   plot(p1l,1:n,'white.',p1r,1:n,'white.','LineWidth',2);
%   plot(p2l,1:n,'white.',p2r,1:n,'white.','LineWidth',2);
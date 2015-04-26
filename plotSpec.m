%function plotSpec(ret,t0,t1,i)
function plotSpec(ret,t0,t1,i)

% i=1;
sr=ret.sr;


x=ret.x(i,:);
t=(0:length(x)-1)/sr;
% figure(3);
% plot(t,x);
% xlim([t0 t1]);

if isfield(ret,'plotcolor');
    plotcolor=ret.plotcolor;
else
    plotcolor='b';
end



N=(t1-t0)*sr+1;
ff=linspace(0,sr-1/(t1-t0),N);

fs=abs(fft(x(t0*sr:t1*sr)));
cs=abs(cmtm2(x(t0*sr:t1*sr),x(t0*sr:t1*sr),2));
[ps w]=pmtm(x(t0*sr:t1*sr),2);
ps=abs(ps);
w=w*sr/2/pi;
% figure(4);
% plot(ff,fs/max(fs),'r',ff,cs/max(cs),'g',w,ps/max(ps),'b');
semilogy(w,ps,plotcolor);
 xlim([0 ret.sr/2]);
xlabel('frequency(Hz)');
ylabel('normalized power');
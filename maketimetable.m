
r=0.1:0.1:180;
t=zeros(1,length(r));
for i=1:length(r)
rr=num2str(r(i))
[status,result]=system(['phtimes 10 ' rr ' P']);
t(i)=str2double(result(11:19));
end
clear rr;
rr=r;
tt=t;
save ptimes180 rr tt;
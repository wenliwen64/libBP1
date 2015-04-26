!ls *SAC > saclist
fid=fopen('saclist');

while ~feof(fid)
      dataname=strtrim(fgetl(fid));
      remain=dataname;
      for i=1:2
     [token, remain] = strtok(remain,'.');
      end
%     s=textscan(fid,'%s',12,'delimiter','.');
     system(['ls *' token '*SAC > ' token '_list' ]);
    disp(['station ' token]);
end
!ls *_list > stanamelist
fid=fopen('stanamelist');
while ~feof(fid)
    dataname=strtrim(fgetl(fid));
%     system(['ls *' dataname '*SAC > tmp1' ]);
    system(['./merge.sh ' dataname]);
    disp(['station' dataname]);
end

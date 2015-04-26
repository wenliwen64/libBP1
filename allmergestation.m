% ! ls *BHN*SAC *BH1*SAC> saclist
% fid=fopen('saclist');
% 
% while ~feof(fid)
%       dataname=strtrim(fgetl(fid));
%       remain=dataname;
%       for i=1:8
%      [token, remain] = strtok(remain,'.');
%       end
% %     s=textscan(fid,'%s',12,'delimiter','.');
% system(['ls' ' *' token '*BHN*SAC' ' *' token '*BH1*SAC' ' > ' token '.list' ]);
% 
%     disp(['station ' 'BHN ' token]);
% end
% !ls *.list > stanamelist
% 
% fid=fopen('stanamelist');
% while ~feof(fid)
%     dataname=strtrim(fgetl(fid));
% %     system(['ls *' dataname '*SAC > tmp1' ]);
%     system(['~/matlab/formalfunc/merge.sh ' dataname]);
%     disp(['station' dataname]);
% end
% ! ls *.r > filelist
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ! ls *BHE*SAC *BH2*SAC> saclist
% fid=fopen('saclist');
% 
% while ~feof(fid)
%       dataname=strtrim(fgetl(fid));
%       remain=dataname;
%       for i=1:8
%      [token, remain] = strtok(remain,'.');
%       end
% %     s=textscan(fid,'%s',12,'delimiter','.');
%      system(['ls' ' *' token '*BHE*SAC' ' *' token '*BH2*SAC' ' > ' token '.list' ]);
%     disp(['station ' 'BHE' token]);
% end
% !ls *.list > stanamelist
% 
% fid=fopen('stanamelist');
% while ~feof(fid)
%     dataname=strtrim(fgetl(fid));
% %     system(['ls *' dataname '*SAC > tmp1' ]);
%     system(['~/matlab/formalfunc/merge.sh ' dataname]);
%     disp(['station' dataname]);
% end
% ! ls *.r > filelist
% ! rm -f list
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ls *BHZ*SAC> saclist
fid=fopen('saclist');

while ~feof(fid)
      dataname=strtrim(fgetl(fid));
      remain=dataname;
      for i=1:2
     [token, remain] = strtok(remain,'.');
      end
%     s=textscan(fid,'%s',12,'delimiter','.');
     system(['ls *' token '*BHZ*SAC > station' token '.list' ]);
    disp(['station ' 'BHZ ' token]);
end
!ls *.list > stanamelist
!awk '{print "~/Programs/bin/merge.sh", $1}' stanamelist | sh
fid=fopen('stanamelist');
while ~feof(fid)
    dataname=strtrim(fgetl(fid));
%     system(['ls *' dataname '*SAC > tmp1' ]);
    system(['~/Programs/bin/merge.sh ' dataname]);
    disp(['station' dataname]);
end
% ! ls *.r > filelist
% ! rm -f *list
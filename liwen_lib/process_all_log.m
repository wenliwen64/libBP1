filelist = importdata('./logfilelist');
for i=1:numel(filelist)
    readlog_function_nopic(char(filelist(i)));
    i
end

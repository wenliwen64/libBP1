% read all sac files and stored into mat
fid=fopen('/Users/lwen/Documents/Seismology_Research/Nepal_2015/nepal_hinet_201404250605_future/jst_corrected/all/sac.dat')
tline = fgets(fid);
count = 0;
while ischar(tline) 
    count = count + 1;
    num_tline = numel(tline);
    tline = tline(1:num_tline-1);
    S(count) = rdsac(tline);
    tline = fgets(fid);
end

save('Hinet_nepal.mat', '-v7.3', 'S');
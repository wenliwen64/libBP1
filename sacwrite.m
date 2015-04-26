function  sacwrite(filename, sac_struct_data)

%     sacwrite('filename', sac_struct_data)
%     Write an ascii format sruct data (read in matlab using
%     sacread, to a binary sacfile. All the headers for the
%     sacfile are included in 'sac_struct_data'
%     fhdr = fread(fid,14*5,'float');
%     ihdr = fread(fid,8*5,'int');
%     chdr = fread(fid,24*8,'char')';
%     chdr = char(chdr);
%     npts = ihdr(10);
%     data = fread(fid, npts,'float');
     

if(iscell(filename) == 1) 
    filename = char(filename);
end
    
fid = fopen(filename,'w');

SAC = sac_struct_data;
npts = SAC.npts;

names = fieldnames(SAC);
Num = length(names);
nf = 70;
ni = 40;
nc = 23;

for i=1:nf
    fwrite(fid, getfield(SAC,char(names(i))),'float');
end
for k=1:ni
    i = nf+k;
    fwrite(fid, getfield(SAC,char(names(i))),'int');
end
for k=1:nc
    i = nf+ni+k;
        fwrite(fid, getfield(SAC,char(names(i))),'char');
end
% to write in the data
fwrite(fid, getfield(SAC,'data',{1:npts}),'float32');
fclose(fid)

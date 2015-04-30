function sh = rdsac_group(directory)
    char(directory)
    list = ls(char(directory));
    list = strread(list, '%s', 'delimiter', '\n');
    for i = 1:numel(list)
        sh(i) = rdsac(list{i});
    end
    
    save('sh.mat', 'sh', '-v7.3');
end
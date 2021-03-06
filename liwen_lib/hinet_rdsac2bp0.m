% convert sac to mat file. 
% usage: hinet_rdsac2bp0('./*sac', 'result.mat')
function hinet_rdsac2bp0(directory, matname)
    ret = dir(directory);
    
    for i = 1:numel(ret)
        sh(i) = rdsac(ret(i).name);
    end
    save(matname, 'sh', '-v7.3');
    clear all;
end
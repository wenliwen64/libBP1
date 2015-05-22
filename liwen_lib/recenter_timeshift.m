for i = 1:71
    ret1.timeshift(:,:,i) = ret1.timeshift(:,:,i) - ret1.timeshift(21,21,i);
end
 ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
            for ii = 1:6; axes(ha(ii)); pcolor(zeros(4,4)); end
%             set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
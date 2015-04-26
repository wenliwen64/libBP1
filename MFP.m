
n [ Y, r, V, DV] = MFP(P,Lambda,c)

%This function solves multiple fares dynamic program assuming all 
%the demands are Poisson. 
%Input:     (P, Lambda, c)
%   1)  P is the fare prices vector in decreasing order.
%   2)  Lambda is the corresponding demand rates vector. 
%   3)  Capacity c.
%   The number of class is determined by P, extra data in Lambda and c
%   will be ignored. This function does minimal validation checks.
%Output:    Y, r, V
%       Y is a row vector of protection levels $y_1,\ldots, y_{n-1}$ 
%       r is the maximum expected revenue $V_n(c)$.
%       V the matrix $V_j(x)$
%       DV the matrix $\Delta V_j(x)$.
%       The  can write $V, \DeltaV$,  and $Y$ to an Excel file called MF.
%       Writing to Excel is slow so that function is commented.
%       Please remove the comments in section Reporting Output to allow the write functionaltiy.
%   Use it as
%     [Y, r, V, DV]= MFP(P,Lambda,c) if you want all the output or as
%     [Y, r] = MFP(P,Lambda,c) if you only want $Y$ and $r$.


% Last edited: Guillermo Gallego February 2010.
%%   Check parameters
n=length(P);
%if (n< 2) || (length(Lambda) < n) || (nargin <3 )
%    disp('Error: insurficient parameters.');
%elseif (c<=0)
%    disp('Error: Invalid parameters.');
%else
%% Prepare demand distribution. 
% Data for stage i (value j) % is stored in (row i, column j+1).
V =   zeros(n, c+1); 
      for stage = 1:n
            Pdf(stage,1) = poisspdf(0,Lambda(stage));
            Cdf(stage,1) = Pdf(stage,1);
       for x = 2:c+1 
           Pdf(stage, x)= poisspdf(x-1, Lambda(stage));
           Cdf(stage,x) = Cdf(stage,x-1) + Pdf(stage,x);
       end
    end
    
    
 %% Stage 1
   stage = 1;
        V(1,1) = 0;
        DV(1,1) = P(1);
        StageY(1,1) = 0;
           for x=1:c
              DV(1,x+1) = P(1)*(1-Cdf(stage,x)); 
              V(1,x+1) = V(1,x) + DV(1,x+1) ; %V(j,x+1) = V_j(x).
           end;  
           for x = 0:c
                if DV(1,x+1) > P(2)
                  t(1) = x+1; %t(1) = y_1+1 at end of for loop
                end
           end
            
  %% Stage 2 to n-1  
        for stage = 2: n-1
        for x = 1:t(stage-1)
            DV(stage,x) = DV(stage-1,x);
            V(stage,x) = V(stage-1,x);   %since $V_j(x) = V_{j-1}(x) on x \leq y_j
        end
        for x = t(stage-1)+1:c+1;
             j = x-t(stage-1);
                temp = 0;
                for i=0:j-1;
                    temp = temp + DV(stage-1,t(stage-1)+j-i)*Pdf(stage,i+1); 
                end 
             DV(stage,x)= temp + P(stage)*(1-Cdf(stage,j));
             V(stage,x) = V(stage,x-1) + DV(stage,x);
        end
            for x=0:c
                if DV(stage,x+1) > P(stage +1)
                    t(stage) = x+1;
                end
            end
        end
        
        
         %% Stage n
        stage=n;
        for x = 1:t(stage-1)
            V(stage,x) = V(stage-1,x);   %since $V_k(x) = V_1(x) on x \leq y_1
            DV(stage,x) = DV(stage-1,x);
        end
        for x = t(stage-1)+1:c+1;
             j = x-t(stage-1);
                temp = 0;
                for i=0:j-1;
                    temp = temp + DV(stage-1,t(stage-1)+j-i)*Pdf(stage,i+1); 
                end 
             DV(stage,x)= temp + P(stage)*(1-Cdf(stage,j));
             V(stage,x) = V(stage,x-1) + DV(stage,x);
        end

         r=V(n,c+1);

%%  Get optimal policy        
    Y(n)=c;
    Y(1)=t(1)-1;
    for stage=2:n-1
        Y(stage)=max( t(stage)-1, Y(stage-1));
    end
   
    B= V;
    DB = B(:,1);
    for x =2:c+1
    DB(:,x) = B(:,x) - B(:,x-1);
    end
    
       
    DV = DB;
    B = B';
    DB = DB';


%Preparing Output
stg=1:1:n;
cap=0:1:c;
cap = cap';

%Reporting Output
%tit1={'V_j(x)'};
%flag1=xlswrite('MFP.xls',tit1,'V_j(x)','A2');
%info1=[B];
%flag2=xlswrite('MFP.xls',info1, 'V_j(x)', 'B3');
%info2=[stg];
%flag3=xlswrite('MFP.xls',info2, 'V_j(x)', 'B2');
%info3=[cap];
%flag3=xlswrite('MFP.xls',info3, 'V_j(x)', 'A3');

%tit1={'DeltaV_j(x)'};
%flag1=xlswrite('MFP.xls',tit1,'DeltaV_j(x)','A2');
%info1=[DB];
%flag2=xlswrite('MFP.xls',info1, 'DeltaV_j(x)', 'B3');
%info2=[stg];
%flag3=xlswrite('MFP.xls',info2, 'DeltaV_j(x)', 'B2');
%info3=[cap];
%flag3=xlswrite('MFP.xls',info3, 'DeltaV_j(x)', 'A3');


%tit1={'y'};
%flag1=xlswrite('MFP.xls',tit1,'Y','A2');
%info1=[Y];
%flag2=xlswrite('MFP.xls',info1, 'y', 'B2');




end



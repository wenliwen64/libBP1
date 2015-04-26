function x1= randphase(x)
    n=length(x);
    tmpf=fft(x);
    tmpa=abs(tmpf);
    tmpph=tmpf./tmpa;
    
  
    if mod(n,2)==1
    r=2*pi*rand(1,(n-1)/2);% making noise uniform distribution
%         r=pi*randn(1,(n-1)/2);% making noise gaussian distribution with sigma=pi

    tmpph(2:(end-1)/2+1)=cos(r)+1i*sin(r);
%         tmpph(2:(end-1)/2+1)=tmpph(2:(end-1)/2+1)+cos(r)+1i*sin(r);

    tmpph((end-1)/2+2:end)=fliplr(tmpph(2:(end-1)/2+1))';
 
    else
    r=2*pi*rand(1,(n-2)/2);% making noise
%             r=pi*randn(1,(n-1)/2);% making noise gaussian distribution with sigma=pi

     tmpph(2:end/2)=cos(r)+1i*sin(r);
%           tmpph(2:end/2)=tmpph(2:end/2)+cos(r)+1i*sin(r);

    tmpph(end/2+2:end)=fliplr(tmpph(2:end/2))';   
    end
    
           x1=ifft(tmpph.*tmpa);

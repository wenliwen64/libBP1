function xf=bandpass(x,Fs,colo,cohi,npol,npas,tipe)
% xf=BANDPASS(x,Fs,colo,cohi,npol,npass,tipe)
%
% Filters signal 'x' with filter 'tipe' and corner
% frequencies 'cohi' and 'cohi' in Hz with 'npol' the
% number of poles and in 'npas' passes. Sampling rate is 'Fs' in Hz.
%
% Compare in SAC bp butter co 0.05 5 n 2 p 1
%
% Last modified by fjsimons@alum.mit.edu, Jan 09, 2003

defval('npol',2)
defval('npas',1)
defval('colo',0.05)
defval('cohi',0.50)
defval('Fs',110)
defval('tipe','butter')

disp(sprintf('BANDPASS %3.3f Hz %i pass %i poles %s',co,npas,npol,tipe))
                                        
% Corner frequency is in Hertz, now it is as a fraction of
% half the sampling rate.
Wn=2*[colo cohi]/Fs;

[B,A]=feval(tipe,npol,Wn);

xf=filter(B,A,detrend(x(:)));

if npas==2
  xf=flipud(filter(B,A,detrend(flipud(xf(:))))); 
end

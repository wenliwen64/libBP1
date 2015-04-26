function xf=hipass(x,Fs,co,npol,npas,tipe)
% xf=HIPASS(x,Fs,co,npol,npass,tipe)
%
% Filters signal 'x' with filter 'tipe' and corner
% frequency 'co' in Hz with 'npol' number of poles and in
% 'npas' passes. Sampling rate is 'Fs' in Hz.
%
% Compare in SAC hp butter co 5 n 2 p 1
% See FilterTest for comparison.
%
% See also LOWPASS, BANDPASS
%
% Last modified by fjsimons@alum.mit.edu, Jan 09, 2003

defval('npol',2)
defval('npas',1)
defval('co',5)
defval('Fs',110)
defval('tipe','butter')

disp(sprintf('HIPASS %3.3f Hz %i pass %i poles %s',co,npas,npol,tipe))
                                        
% Corner frequency is in Hertz, now it is as a fraction of
% half the sampling rate.
Wn=2*co/Fs;

[B,A]=feval(tipe,npol,Wn,'high');

xf=filter(B,A,detrend(x(:)));

if npas==2
  xf=flipud(filter(B,A,detrend(flipud(xf(:))))); 
end

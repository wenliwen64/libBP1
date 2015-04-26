% XCORRSHIFT Cross-correlation time-delay measurement
%
% SYNTAX	[dsamp,maxc,xcor] = xcorrshift(s1,s2)
% 		[dsamp,maxc,xcor] = xcorrshift(s1,s2,win)
%
% PURPOSE	Slides signal s1 through signal s2 to find
%		the relative time-delay that maximizes the cross-correlation.
%
% INPUT		s1	reference signal (short window) 
%		s2	target signal (long window)    
%		win	window for search of max coherency (in samples) [] 
%
% OUTPUT	dsamp	relative delay (in samples), 
%			sub-sample precision estimated by quadratic interpolation
%			ALWAYS >=0
%			0 MEANS ALIGNED
%		maxc	optimal coherency (can be negative) [-1:1]
%		xcor 	cross-correlation time series = Sxy/sqrt(Sxx*Syy)
%
% NOTE		XCORRSHIFT is usually applied for cross-correlation picking.
%		If s1 and s2 are selected windows of two signals,
%		starting at samples p1 and p2 respectively (raw time-picks),
%		i.e. s1 = sis1(p1:...) and s2 = sis2(p2:...),
%		then the new pick for s2, based on s1, is p2+dsamp.
%		THIS FUNCTION GIVES ONLY POSITIVE DELAYS (dsamp).
%		The begining of window s2 must be such that the signal inside
%		is DELAYED with respect to the reference windowed signal s1:
%		you must leave enough header samples in s2.
%
% NOTE		Cross-correlation is computed by FFT, so XCORRSHIFT is optimal
%		when signal lengths are of the same order.

% Jean-Paul Ampuero	ampuero@erdw.ethz.ch

function [dsamp,maxc,xcor] = xcorrshift(s1,s2,win,only_positive)

n1=length(s1);
n2=length(s2);
if n1>n2, error('Reference signal S1 cannot be longer than target S2'), end
nx = n2-n1+1; % length of non-spurious part of xcor (no wraparound) 
if nargin<3, win=[]; end
if isempty(win)
  win=[1:nx]; 
else
  win=[max(1,win(1)):min(nx,win(end))];
end
if nargin<4, only_positive=0; end

% cross-correlation
nfft = 2^nextpow2(n2+n1);
f1=fft(s1,nfft);
f2=fft(s2,nfft);
xcor=real(ifft( conj(f1).*f2 ));
xcor=xcor(1:nx);

% scale by sqrt(norm(s1)*norm(s2win)) where s2win is the moving window of s2
s2s2 = s2.*s2;
scal=zeros(nx,1);
scal(1) = sum(s2s2(1:n1));
for k=1:nx-1,
  scal(k+1) = scal(k)+ s2s2(n1+k)-s2s2(k);
end
scal = sqrt(scal)*norm(s1);
xcor=xcor./scal ;

% optimal lag index (=delay+1)
if only_positive
  [maxc,maxi]=max(xcor(win));
else
  [maxc,maxi]=max(abs(xcor(win)));
end
maxi = maxi +win(1)-1;
maxc=xcor(maxi);

% sub-sample precision
if maxi>1 & maxi<nx-1
  xc1 = 0.5*( xcor(maxi+1)-xcor(maxi-1) );
  xc2 = xcor(maxi-1)-2*xcor(maxi)+xcor(maxi+1);
  maxi = maxi-xc1/xc2;
  maxc = maxc-0.5*xc1*xc1/xc2;
end

% lag --> delay
dsamp= maxi-1;

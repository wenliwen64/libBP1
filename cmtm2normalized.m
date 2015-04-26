%[Coh,Pha,DPha] = cmtm(X,Y,NW,K,freq)
%
%CMTM coherence function estimate 
%via the Thomson multitaper method (MTM).
%
% INPUT                X,Y        equal length (2^p) signals
%                NW        time-bandwidth product
%
% OUTPUT        Coh        double sided coherence spectrum
%                Pha        phase cross-spectrum
%                DPha        variance of phase cross-spectrum
%
% See also PMTM and Frederik Simons' PMTM2D
%
function [Coh,Pha,DPha] = cmtm2(X,Y,NW,K,freq,FIG)

persistent DPSS_K DPSS_NT DPSS_NW DPSS_E DPSS_V

if ~exist('FIG','var'), FIG=0; end
if ~exist('K','var'), K=[]; end
if ~exist('freq','var'), freq=[1:length(X)]; end

X=X(:)-mean(X);        % make input column vectors and demean 
Y=Y(:)-mean(Y);
nt = length(X);
nf = length(freq);
if isempty(K), K=max(round(2*NW)-1,1); end

% Compute the data windows (columns of E)
if isempty(DPSS_V) | nt~=DPSS_NT | NW~=DPSS_NW | K~=DPSS_K
  disp('Generating DPSS data')
  DPSS_NW=NW;
  DPSS_NT=nt;
  DPSS_K=max(round(2*DPSS_NW)-1,1); % Shannon number-1
  DPSS_K=min(DPSS_K,K);
  [DPSS_E,DPSS_V]=dpss(nt,NW,DPSS_K);
  DPSS_E=DPSS_E(:,1:DPSS_K); 
  DPSS_V=DPSS_V(1:DPSS_K);
  disp(sprintf('DPSS: NW=%u  K=%u  NT=%u',DPSS_NW,DPSS_K,DPSS_NT))
end

% Spectral and cross-spectral densities
Xwigs=fft(DPSS_E.*repmat(X,1,DPSS_K),nt);
Xwigs=Xwigs(freq,:);
SXX=(Xwigs.*conj(Xwigs))*DPSS_V;
  
Ywigs=fft(DPSS_E.*repmat(Y,1,DPSS_K),nt);
Ywigs=Ywigs(freq,:);
SYY=(Ywigs.*conj(Ywigs))*DPSS_V;

SXY=(Xwigs.*conj(Ywigs))*DPSS_V;
    
% Coherence spectrum
Coh = SXY;%./ sqrt((SXX.*SYY));
%Coh = (SXY.*conj(SXY)) ./ (SXX.*SYY) ;
%Coh = sqrt(Coh);

% Phase cross-spectrum
if nargout>1
  Pha = atan2(imag(SXY), real(SXY));
  if FIG
    oldfig=gcf;
    figure(FIG)
    plot(freq,Pha,'k-')
    hold on
  end
  %jackknife variance estimate for the phase spectrum
  if nargout>2 & DPSS_K>1 
    DPha=zeros(nf,1);
    for k=1:DPSS_K
      jackv=DPSS_V;
      jackv(k)=0;
      jackv=jackv/sum(jackv);
      SXY=(Xwigs.*conj(Ywigs))*jackv;
      DPha=DPha+SXY./abs(SXY); % sum jk phase factors
      if FIG, plot(freq,atan2(imag(SXY), real(SXY))); end
    end
    
    %Thomson and Chave 1991 p.91 eq.2.63
    DPha=sqrt(2*(DPSS_K-1)*(1-abs(DPha)/DPSS_K));
  else
    DPha=ones(nf,1)*0.1*pi;
  end
  if FIG
    hold off
    figure(oldfig)
  end
  
end

function [Y,spar,Sps] = SpecSyn2(N,samp,corr,acf,pow2,Rseed);
%
%  [Y,spar,Sps] = SpecSyn2(N,samp,corr,'acf',pow2,Rseed) 
%  generates a 2D-random field Y of size (Nz+1 x Nx+1) with
%  possibly variable spatial sampling in z,x-direction. 
%  This function simulates anisotropic random fields, i.e. 
%  the correlation length in both directions can be different; 
%  rectangular grid dimensions are also possible.
%  corr contains the corr. length ax, az and the Hurst exponent H 
%  or the fractal dimension D; the autocorrelation function to 
%  be used has to be specified by acf;
%  NOTE: D = 3-H, and larger D (smaller H) yield "rougher" fields
%
%  The algorithm is based on the paper by 
%
%        Pardo-Iguzquiza, E. and Chica-Olma, M. (1993)
%        The Fourier integral method: and efficient spectral method for 
%        simulation of random fields, Mathematical Geology, 25, p177-217.
%
%  but extends their method to handle rectangular grids and variable
%  sampling in the two directions.
%
%  INPUT:
%  N -            grid dimensions [Nz Nx]
%  samp         - desired sampling, [dz dx]; if only dx, then dz = dx
%  corr  - corr = [az ax]   for 'gs' and 'ex' 
%           corr = [az ax H] for 'ak'; note that 0 <= H <= 1
%           corr = [D kc]    for 'fr'; D is the fractal dimension,
%                  kc the corner wavenumber beyond which the spectrum decays
%  acf   - autocorrelation function: 
%          'gs' or 'GS' - Gaussian
%          'ex' or 'EX' - Exponential
%          'ak' or 'AK' - anisotropic vonKarman
%          'fr' or 'FR' - fractal distribution
%  pow2  - 'y' if FFT on a grid of size of power-of-two in both
%           directions; this requires an adjustment of the sampling
%           interval, the new sampling is saved in the structure samp
%  Rseed - seeds for random number generators; if omitted or empty
%           Rseed = sum(100*clock) is used and returned with output
%           structure spar
%           [Rpseed Rsseed] for the phase and small random spectral part
% 
%  OUTPUT:
%  Y    - real random field whose size is determined by the sample
%          spacing in each direction, the number of points and whether
%          the pow2-option is given or not. The orientation is 
%          Y = [z-dir x-dir] 
%  spar        - structure with length vectors, sampling in both direction
%          and other relevant parameters; it also contains the random
%          seed number, useful to reproduce realizations
%  Sps  - structure containing the power spectrum ("upper half only")
%     and the wavenumber vectors.
%
%  This function is a modification and - hopefully - improvement of a
%  function provided by the SRB-toolbox. It uses the spectral
%  synthesis method, originally coded up by Philippe Rio. 
%
%  Written by Martin Mai (mmai@pangea.Stanford.EDU) 
%  07/16/98
%  last change 03/01/2000
% ------------------------------------------------

check = 'n';                %% set to 'y' if you want to create
                        %% a simple out put plot to check the
                        %% the spectra and the resulting field

%%% check input variables
if nargin < 4; error('Not enough input arguments');
  elseif nargin == 4; pow2 = 'n'; Rseed = [];
  elseif nargin == 5, Rseed = []; 
end;
if length(samp) == 2; Dz = samp(1); Dx = samp(2);
  elseif length(samp) == 1; Dz = samp(1); Dx = Dz; 
end
if length(N) == 1, N = [N N]; end

%%% set size for spectral synthesis tool that generates
%%% fields of size (2*rmz+1, 2*rmx+1);
Nx = N(2); Nz = N(1);
px = ceil(Nx/Dx); pz = ceil(Nz/Dz);
if mod(Nx,2) == 1; nx  = (Nx+1)/2; else nx  = Nx/2; end
if mod(Nz,2) == 1; nz  = (Nz+1)/2; else nz  = Nz/2; end
if mod(px,2) == 1, rmx = (px+1)/2; else rmx = px/2; end
if mod(pz,2) == 1, rmz = (pz+1)/2; else rmz = pz/2; end

if (Dx==1) & (mod(Nx,2) == 1); rmx = (Nx-1)/2; end
if (Dz==1) & (mod(Nz,2) == 1); rmz = (Nz-1)/2; end


%%% adjust sampling, if necessary
if (2*rmx+1)*Dx ~= Nx+Dx
        dx=Nx/(2*rmx);
        xs = sprintf('%.3f --> %.3f',Dx,dx);
        disp(['  ** need to adjust x-sampling: ',xs]);
else
        dx=Dx;
end
if (2*rmz+1)*Dz ~= Nz+Dz
        dz=Nz/(2*rmz);
        zs = sprintf('%.3f --> %.3f',Dz,dz);
        disp(['  ** need to adjust z-sampling: ',zs]);
else
        dz=Dz;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% this check didn't work quite okay, but I'll leave it in
%%%% in case I need to refer back to it later
%%if mod(Nx,2) == 1; nx = (Nx+1)/2; else nx = Nx/2; end
%%if mod(Nz,2) == 1; nz = (Nz+1)/2; else nz = Nz/2; end
%%rmx = nx/Dx; rmx2 = round(rmx); dx = nx/rmx2;
%%rmz = nz/Dz; rmz2 = round(rmz); dz = nz/rmz2;
%%% print adjusted sampling to screen
%%if (rmx-rmx2) ~= 0;
%%        rmx = rmx2; 
%%        xs = sprintf('%.3f --> %.3f',Dx,dx);
%%        disp(['  ** need to adjust x-sampling: ',xs]);
%%end
%%if (rmz-rmz2) ~= 0;
%%        rmz = rmz2;
%%        zs = sprintf('%.3f --> %.3f',Dz,dz);
%%        disp(['  ** need to adjust z-sampling: ',zs]);        
%%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% now, the resulting field would be of size (2*rmz+1,2*rmx+1)
%%% which in most cases won't be close to a power of two;
%%% to speed up the FFT, we can force it to be a power of two
%%% which requires us to redefine the sample spacing
if pow2 == 'y'
        px = nextpow2(4*rmx);
        pz = nextpow2(4*rmz);
        dx = Nx/(2^px);
        dz = Nz/(2^pz);
        rmx = ceil(nx/dx);
        rmz = ceil(nz/dz);
end

%%% wavenumber vector in [-pi,pi]
kx = (1/dx)*(2*pi)*(1/(2*nx+1))*[-nx:dx:nx];
kz = (1/dz)*(2*pi)*(1/(2*nz+1))*[-nz:dz:nz];

%%% get some variables from input variable corr
if  (acf == 'fr') | (acf == 'FR')
   if (length(corr) == 2)
        D = corr(1); kc = corr(2);
   else
        D = corr(1); kc = 0; 
        disp('  ** Corner wavenumber kc set to zero **')
   end
elseif (acf == 'ex') | (acf == 'EX') | (acf == 'gs' ) | (acf == 'GS')
        ax = corr(2);  az = corr(1);
elseif (acf == 'ak') | (acf == 'AK')
        ax = corr(2);  az = corr(1); H = corr(3);
end

%%% wavenumber field, for two of the four quadrants needed in the IFFT
kr = zeros(rmz+1,2*rmx+1);
for j = 1:(2*rmx+1)
    for i = 1:(rmz+1)
        if ((acf == 'fr') | (acf == 'FR'))
                kr(i,j) = sqrt((kz(i)^2) + (kx(j)^2));
        else
                  kr(i,j) = sqrt((az^2)*(kz(i)^2) + (ax^2)*(kx(j)^2));
        end
    end
end

%%% calculate power spectral density, depending on selected ACF
if ((acf == 'gs') | (acf == 'GS')) 
        PS = 0.5*ax*az*exp(-0.25*kr.^2);
elseif ((acf == 'ex') | (acf == 'EX'))
          PS = (ax * az)./(1 + (kr.^2)).^1.5;
elseif ((acf == 'ak') | (acf == 'AK'))
        k = kr(:); k = k(k>0);
        coef = 4*pi*H*ax*az./besselk(H,min(k));
           PS = coef./(1 + (kr.^2)).^(H+1);
elseif ((acf == 'fr') | (acf == 'FR'))
        decay = 0.5*(8-2*D);
         %% to ensure proper scaling of the power spectrum 
        %% in k-space we cannot allow kr == 0
        if min(kr(:)) == 0
          [p,q] = find(kr == 0);
          k(p,q) = mean(mean(kr(p-1:p,q-1:q)));
        end
        %% take care of corner wavenumber                        
        kr(kr <= pi*kc) = pi*kc;                
        PS = 1./((kr.^2).^decay);
        %%kxpos = kx(kx>=0);
        %%ix = find(kxpos <= pi*kc);
        %%kcval = 1/(((kxpos(ix(end)))^2)^decay);
               %%PS(kr <= pi*kc) = kcval;        
end

%%% the IFFT needs the spectrum normalized to max.amp unity
PS = PS./max(PS(:));
AM = sqrt(PS);

%%% --- now deal with the random phase ---
%%% initialize random number generator
if isempty(Rseed) == 1;
  Rseed = zeros(1,2);
  Rseed(1) = sum(100*clock);
  Rseed(2) = sum(109*clock);
elseif length(Rseed) == 1
    Rseed = [Rseed 2*Rseed];
end
rand('seed',Rseed(1));
randn('seed',Rseed(2));

%%% random phase in [0,pi]
PH = 2*pi*rand(size(kr));         

%%% assemble the random field in FT-domain        
% add small random high-wavenumber components
%x = (1 + 0.5*randn(size(kr)));
x = 1;
RAD = AM.*x;

%%% set DC-component to different value, if desired
%%% NOTE that this only changes the overall 'level' of
%%% the field, not its appearance, but has significant
%%% impact on the Fourier transform which reflects the
%%% DC-value at the smallest wavenumber
Neff = 0;
RAD(rmz+1,rmx+1) = Neff;                % "nugget" effect, zero-mean field
AM(rmz+1,rmx+1)  = RAD(rmz+1,rmx+1);
Y = RAD.*cos(PH) + sqrt(-1)*RAD.*sin(PH);

%%% the following takes care of the conjugate symmetry condition
%%% in order to ensure the final random field is purely real
U = zeros(2*rmz+1,2*rmx+1);        % becomes the conjugate symmetric field
Y = [Y;conj(fliplr(flipud(Y(1:rmz,:))))];
for i = 1:rmx
     Y(rmz+1,-i+2*rmx+2) = conj(Y(rmz+1,i));
end
for i = 1:rmz+1
   for j = 1:rmx+1
      U(i,j) = Y(i+rmz,j+rmx);
    end
end
for i = rmz+2:2*rmz+1
   for j = rmx+2:2*rmx+1
      U(i,j) = Y(i-rmz-1,j-rmx-1);
    end
end
for i = 1:rmz+1
   for j = rmx+2:2*rmx+1
      U(i,j) = Y(i+rmz,j-rmx-1);
    end
end
for i = rmz+2:2*rmz+1
   for j = 1:rmx+1
      U(i,j) = Y(i-1-rmz,j+rmx);
    end
end

%%% take 2D-inverse FFT to obtain spatial field; imaginary parts of order 1e-13
%%% due to machine precision are removed; remove mean and scale to unit variance
if pow2 == 'y'
        Y = real(ifft2(U,4*rmz,4*rmx));
else
        Y = real(ifft2(U));
end
Y = Y/std2(Y);                                % standard deviation of unity
%if mean(Y(:)) < 0; Y = (-1)*Y; end        % have mean being positive (though small)
                                        
%%% final output structure with length vectors,sample spacing
%%% and model parameters
spar.dim = N;
spar.samp = [Dz Dx];
spar.sampM = [dz dx];
spar.lx = [0:dx:dx*(size(Y,2)-1)]; 
spar.lz = [0:dz:dz*(size(Y,1)-1)]; 
spar.size = size(Y);
spar.corr = corr;
spar.acf = acf;
spar.pow2 = pow2;
spar.Rseed = Rseed;

%% assemble structure with spectra in pos. quadrant
kpx = 1:size(PS,2);
kpz = 1:size(PS,1);
Sps.PD  = PS(:,kpx);
Sps.kpx = kx(kpx)';
Sps.kpz = kz(kpz)';
Sps.PDx = Sps.PD(end,:)';
Sps.PDz = Sps.PD(:,1);

%%% plot spectrum for error checking
if check == 'y';

figure
fullpage
subplot(221);
imagesc(Sps.kpx,Sps.kpz,log10(Sps.PD)); colorbar
xlabel('kx'); ylabel('kz');
title('2D-Power Spectral Density (in log-units)','FontS',12);
axis equal; axis tight;


subplot(222)
loglog(-Sps.kpz,Sps.PDz,'r','LineW',2); hold on
loglog(Sps.kpx,Sps.PDx,'b','LineW',2);
if acf == 'fr'
   loglog(pi*kc*ones(1,10),linspace(1,min(Sps.PD(:)),10),'k--');
   tstr = ['kc = ',num2str(kc),', dz, dx = ',...
        num2str(samp(1)),', ',num2str(samp(2))];
else
  tstr = ['az,ax = ',num2str(az),', ',num2str(ax),', dz, dx = ',...
        num2str(samp(1)),', ',num2str(samp(2))];
end
legend('z','x',3);
xlabel('wavenumber');
ylabel('psd');
title(['1D-spectra (',tstr,')'],'FontS',12);
axis([0 max([max(Sps.kpx) max(-Sps.kpz)]) min(Sps.PD(:)) 2]);
axis square
hold off

subplot(212);
imagesc(spar.lx,spar.lz,Y); colorbar
xlabel('X'); ylabel('X');
title('Resulting Random Field','FontS',12,'FontW','bo');
axis equal; axis tight

end

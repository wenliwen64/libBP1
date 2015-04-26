function x1=hilshift(x,theta)%phase shift the profile with theta

x1=real(hilbert(x)*exp(-1i*theta/180*pi));

 
close all;
clear all;
lengthx=100;
lengthy=20;
N=[lengthx lengthy];
samp=4;
corr=[1 1 0.8];
acf='ak';
pow2='y';
scalex=128;%lengthx/samp*4;
scaley=32;%lengthy/samp*4;

[Y,spar,Sps]=SpecSyn2(N,samp,corr,acf,pow2);
for i=1:scalex
    lon0(i,:)=linspace(0,1,scaley)*lengthy;
end
for i=1:scaley
    lat0(:,i)=linspace(0,1,scalex)*lengthx;
end
figure(1);
pcolor(lon0,lat0,Y);
colorbar;
shading interp;

P=[7.29320740145629e-60,-2.25713355056818e-54,3.15780784770224e-49,-2.64134885422995e-44,1.47242246990781e-39,-5.77410811052121e-35,1.63959228354946e-30,-3.42090336118034e-26,5.27212978329383e-22,-5.99116415754033e-18,4.98013313198646e-14,-2.98197497509444e-10,1.25141562627232e-06,-0.00349386067308002,5.89649663355069,1273.93757459125;]

%parkersfield velocity model

CFUN= @(b) P(16)+P(15)*b+P(14)*b.^2+P(13)*b.^3+P(12)*b.^4+P(11)*b.^5+(P(10)+P(9)*b+P(8)*b.^2+P(7)*b.^3+P(6)*b.^4+P(5)*b.^5+P(4)*b.^6+P(3)*b.^7+P(2)*b.^8+P(1)*b.^9).*b.^6;
z=linspace(0,20,scalex)*1e3;
for i=1:scalex
    vel(i,:)=ones(1,scaley)*CFUN(z(i));
end
Yv=(1+0.05*Y).*vel;
kk=1;
for i=1:scalex
    for j=1:scaley
        velp(kk)=Yv(i,j);
        kk=kk+1;
    end
end
% generates a heterogeneous rho, vp ,vs model in a regular grid
% to be used in the DIST_HETE1 blocks

ncol=3;
nx=scalex;
nz=scaley;
x0=-50000;
z0=-10000;
dx=lengthx/scalex*1e3;
dz=lengthy/scaley*1e3;
 
% heterogeneous material model: background + random
rho = 2400*ones(1,nx*nz);
vp  = velp;
vs  = 2000/3800*vp;

fid=fopen('material.tab','w');
fprintf(fid,'%u %u %u %f %f %f %f\n',ncol,nx,nz,x0,z0,dx,dz);
fprintf(fid,'%f %f %f\n',[rho;vp;vs]);
fclose(fid);
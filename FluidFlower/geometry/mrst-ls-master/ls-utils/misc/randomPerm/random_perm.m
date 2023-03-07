function [perm,var_lnk_actual]= random_perm(var_lnk,corr_lenx,corr_leny,Nx,Ny,Lx,Ly)

% Wave numbers
kx = (2*pi/Lx)*[0:(Nx/2-1) (-Nx/2):(-1)]';
ky = (2*pi/Ly)*[0:(Ny/2-1) (-Ny/2):(-1)]';
[KX,KY]= meshgrid(kx,ky);
KX2= KX.^2;
KY2= KY.^2;
dkx=sqrt(KX2(1,2)); dky=sqrt(KY2(2,1));


%%%%Spectral density functions%%%%

%Whittle-A spectrum, isotropic medium
%a=pi/(4*sqrt(corr_lenx*corr_leny));
%S=(2*var_lnk*a*a/pi)*( (KX2+KY2) ./ ((a*a+KX2+KY2).^3) );  

%Gaussian spectrum, Random Field Generator, Ruan and McLaughlin, Adv Water Res. 1998,
%S=(1/2/pi)*var_lnk*corr_lenx*corr_leny*exp(-0.5*(corr_lenx^2)*(KX2) - 0.5*(corr_leny^2)*(KY2));

%Spectrum of modified exp autocovariance, Gelhar & Axness, WRR 1983
S=(2/pi)*var_lnk*corr_lenx*corr_leny;
S=S./( (1+(corr_lenx^2).*KX2+(corr_leny^2).*KY2).^3 ); 


H=S.^(0.5); 
theta=2*pi*rand(Ny,Nx); %random phase angle

dZ=H.*exp(1i*theta).*sqrt(dkx*dky)*Nx*Ny;
lnK=ifft2(dZ); %log(perm)
K1=sqrt(2)*real(lnK); %K2=sqrt(2)*imag(lnK); 
var_lnk_actual=var(reshape(K1,Ny*Nx,1));

perm=exp(K1); 


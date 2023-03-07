clc;
%clear all; 
close all;

% Grid
Nx = 1000;
Ny = Nx/2;
Lx = 1;
Ly = 0.5;
dx = Lx/Nx;
dy = Ly/Ny;
x = dx/2:dx:Lx-dx/2;
y = dy/2:dy:Ly-dy/2;
[xx,yy] = meshgrid(x,y);

% Parameters
kmean = 30;
var_lnk = log(40*2);
corr_lenx= 0.15*Lx;
corr_leny= 0.05*Ly;

for n=1:100
[perm,var_lnk_actual]= random_perm(var_lnk,corr_lenx,corr_leny,Nx,Ny,Lx,Ly);
perm = perm + (kmean - mean(perm(:)));
sigma(n) = std(perm(:));
minv(n) = min(perm(:));
maxv(n) = max(perm(:));
meanv(n) = mean(perm(:));
end

surf(xx,yy,perm,'facecolor','interp','edgecolor','none','facelighting','phong');
view([0,0,1]); axis equal tight; 
colorbar; colormap(turbo); %caxis([0 1]); 

disp(['std: ' num2str(mean(sigma))])
disp(['mean: ' num2str(mean(meanv))])
disp(['min: ' num2str(mean(minv))])
disp(['max: ' num2str(mean(maxv))])
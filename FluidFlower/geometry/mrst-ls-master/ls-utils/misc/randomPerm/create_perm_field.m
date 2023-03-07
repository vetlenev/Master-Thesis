clc;clear all; close all;

% Grid
Nx = 10;
Ny = 150;
Lx = 2;
Ly = 150;
dx = Lx/Nx;
dy = Ly/Ny;
x = dx/2:dx:Lx-dx/2;
y = dy/2:dy:Ly-dy/2;
[xx,yy] = meshgrid(x,y);

% Parameters
var_lnk = 0.21;
corr_lenx= Lx;
corr_leny= 0.5*Ly;

[perm,var_lnk_actual]= myrandom_perm(var_lnk,corr_lenx,corr_leny,Nx,Ny,Lx,Ly);

surf(xx,yy,log10(perm),'facecolor','interp','edgecolor','none','facelighting','phong');
view([0,0,1]); %axis equal tight; 
colorbar; colormap(jet); %caxis([0 1]); 

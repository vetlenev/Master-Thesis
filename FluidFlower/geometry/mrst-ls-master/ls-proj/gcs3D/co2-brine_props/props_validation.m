%% Validation of pvtBrineWithCO2BlackOil
% Here, we validate our implementation by comparing the output with
% published results. Specifically, we compare our output with:
%   Fig. 4 and Fig. 5 in Hassanzadeh et al, IGGC (2008).
%   Fig. 2 in Spycher and Pruess, GGA (2005)
%    
% NOTE: Although Fig. 5 in Hassanzadeh et al. (2008) states 50C, the plots 
% correspond to T=45C (Hassan Hassanzadeh, personal communication, 2019).
%
% Documentation on the functionality itself is provided within
% pvtBrineWithCO2BlackOil.
%
clc, clear, close all

% Function inputs
T = {'C', 45};                      % temperature
P = {'MPa','mMn',[0.1, 50, 100]};   % pressure range
S = {'ppm', 'brine', 1.5*10^5};     % salinity
saltVar = true;                     % see pvtBrineWithCO2BlackOil
iteration = false;                  % "
figs = true;                        % set to false if you don't want figs

% Generate properties and plot figures.
% Compare Fig. 1 and Fig. 4 with fig 5c and Fig. 5a in Hassanzadeh et al.
% (2008), Fig. 2 with Fig.2 in Spycher and Pruess (2005), Fig. 3 with 
% Fig. 4 in Hassanzadeh et al. (2008), and Fig. 5 with Fig. 5d in 
% Hassanzadeh et al. (2008) (T=45 C)
t = pvtBrineWithCO2BlackOil(T, P, S, saltVar, iteration, figs);
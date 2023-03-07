function mu_b_co2 = viscBrineCO2IC2012(T, P, m_nacl, w_co2)
%
% DESCRIPTION:
% Calculate the dynamic viscosity of a solution of H2O + NaCl (brine) with
% dissolved CO2.
%
% RANGE:
% The model of Mao & Duan (2009) for brine viscosity reaches 623K, 1000 bar
% and high ionic strength. However, the model used to determine the viscosity 
% when co2 dissolves in the brine (Islam & Carlson, 2012) is based on
% experimental data by Bando et al. (2004) and Fleury and Deschamps (2008),
% who provided experimental data up to P = 200 bar, T extrapolated to
% 100 C, and maximum salinity of 2.7M.
%
% INPUT:
% T:          Double with temperature value in Kelvin
% P:          Double with pressure value in bar
% m_nacl:     Salt molality (NaCl) in mol/kg solvent
% w_co2:      Mass fraction of CO2 in the aqueous solution (i.e. brine)
%
% OUTPUT: 
% mu_b_co2:   Double with dynamic viscosity in Pa*s
%

% Check range
if P > 200
    warning(['Pressure out of measured range in ', mfilename])
end
if T < 273.15 + 35 || T > 373.15
    warning(['Temperature out of tested range in ', mfilename])
end
if m_nacl > 3.1 % m approx, limit of experimental data is 2.738M
    warning(['Salinity of tested range in ', mfilename])
end

% Units
P_mpa = P/10;

% Pure water density (see Islam & Carlson, 2012)
a = 1.34136579*10^2;
b = [-4.07743800*10^3, 1.63192756*10^4, 1.37091355*10^3];
c = [-5.56126409*10^-3, -1.07149234*10^-2, -5.46294495*10^-4];
d = [4.45861703*10^-1, -4.51029739*10^-4];
cp = [1, 2];

rho_h2o = a + sum(b.*10.^(c.*T)) + sum(d.*P_mpa.^cp);                     % [kg/m3]
rho_h2o = rho_h2o/10^3;                                                   % [g/cm^3]

% Pure water viscosity (Mao & Duan, 2009)
d = [0.28853170*10^7, -0.11072577*10^5, -0.90834095*10, ...
     0.30925651*10^-1, -0.27407100*10^-4, -0.19283851*10^7, ...
     0.56216046*10^4, 0.13827250*10^2, -0.47609523*10^-1, ...
     0.35545041*10^-4];
coefs1 = (1:5)-3;
coefs2 = (6:10)-8;

mu_h2o = exp(sum(d(1:5).*T.^coefs1) + sum(d(6:10).*rho_h2o.*T.^coefs2));

% Brine viscosity (H2O + NaCl) ( Mao & Duan, 2009)
A = -0.21319213 + 0.13651589*10^(-2)*T -0.12191756*10^(-5)*T^2;
B = 0.69161945*10^-1 - 0.27292263*10^(-3)*T + 0.20852448*10^(-6)*T^2;
C = -0.25988855*10^-2 + 0.77989227*10^(-5)*T;

mu_b = exp(A*m_nacl + B*m_nacl^2 + C*m_nacl^3) * mu_h2o;

% Brine viscosity with dissolved CO2 (H2O + NaCl + CO2) (Islam & Carlson,
% 2012)
mu_b_co2 = mu_b*(1 + 4.65*w_co2^1.0134);

end
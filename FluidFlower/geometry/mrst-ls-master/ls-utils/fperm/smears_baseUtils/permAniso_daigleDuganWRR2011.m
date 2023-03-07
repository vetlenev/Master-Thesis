%% Permeability anisotropy of clay smears
% Model based on Daigle & Dugan, WRR (2011)
% Assumptions:
%
%   - phi (porosity at the beginning of shear) can be computed from a
%     compaction curve and zf (faulting depth). 
%
%   - theta0: the grain orientation angle in degrees (0 = completely 
%     horizontal, 45 is the initial value right after deposition if they 
%     are not oriented) at the end of compaction and before shearing. It 
%     should in theory be computed by using the uniaxial strain \epsilon_v. 
%     If we assume that most of the soil compressibility comes from the 
%     pore compressibility, we can then calculate \epsilon_v for the soil
%     (estimating the depositional porosity) and use that to compute 
%     theta0. 
%     Alternatively, we can consider theta0 unknown. Since theta (the final
%     grain orientation after shearing) mainly depends on the shear strain
%     (gamma) when gamma (= smear.L/smear.T) > 5, it does not really matter. 
%     So, we pick a value of 30 which corresponds to moderate orientation, 
%     and is likely above that of sediments faulted at shallow depths.
%
%

clear

f       = [0, 1];   % fraction of 
m       = [10, 10];
t       = 100;
T       = 10;
theta0  = 44;
phi     = 0.5;

gamma   = t/T;
theta   = acotd(gamma + cotd(theta0));
meq     = sqrt(1/(f(1)/m(1)^2 + f(2)/m(2)^2));
num     = 1 + ((8*meq/9)*cosd(theta) + (2/pi)*sind(theta))/(3*pi/(8*(1-phi)) - 0.5); 
den     = 1 + ((8*meq/9)*sind(theta) + (2/pi)*cosd(theta))/(3*pi/(8*(1-phi)) - 0.5);
kprime  = (num/den)^2;
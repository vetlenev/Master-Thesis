function [a_M, a_R, sG] = getGasSatFromHeightMob2(T, t, B, b, h, h_T, h_B, swr, snr, varargin)
    z_RT = t + h_T; % bottom of residual part on TOP of col
    z_RB = t + h_B; % top of residual part on BOTTOM of col
    z_M = t + h; % mobile height
    H = B - T;
    
    a_M = max(min((z_M - T)./H, 1), 0); % scaling factor for plume saturation
    a_RB = max(min((B - z_RB)./H, 1), 0); % changed from B - h_B to B - z_RB
    a_RT = 1 - a_M - min(max((B - z_RT)./H, 0), 1); % scaling factor for residual saturation
    
    a_R = max(a_RT + a_RB, 0).*(snr > 0); % max to avoid potential negative saturations
    sG = a_M.*(1 - swr) + a_R.*snr; % if below residual CO2 zone, sG = 0   
    
    if nargin > 9
        % Specific reconstruction for NVE cells
        cNVE = varargin{1};              
        z_RT = b(cNVE) - h_max(cNVE); % fill from bottom
        a_M(cNVE) = 0; % no mobile plume
        a_R(cNVE) = min(max((B(cNVE) - z_RT)./H(cNVE), 0), 1);        
        sG(cNVE) = a_R(cNVE).*snr;
        
        % Specific reconstruction for NVE mobile cells
        cNVEMob = varargin{2}; 
        z_RT = b(cNVEMob) - (h_max(cNVEMob) - h_p(cNVEMob)); % b - residual_at_bottom 
        % z_M same as before (mobile plume filled from top (VE assumption)
        a_R(cNVEMob) = min(max((B(cNVEMob) - z_RT)./H(cNVEMob), 0), 1);
        sG(cNVEMob) = a_M(cNVEMob).*(1-swr) + a_R(cNVEMob).*snr;
    end
end
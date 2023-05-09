function [a_M, a_R, sG] = getGasSatFromHeightVirtual(T, t, B, b, h, h_T, h_B, swr, snr, varargin)
    % Evaluate gas saturation in virtual cell of RVE zone.
    % INPUTS:
    %   T: top depth of virtual cell
    %   t: top depth of parent cell
    %   B: bottom depth of virtual cell
    %   b: bottom depth of parent cell
    %   h: height of mobile plume from parent cell
    %   h_T: max height reached for top segregated plume
    %   h_B: net height reached for plumes from sealing layers
    %   swr: residual water saturation
    %   snr: residual gas saturation
    % RETURNS:
    %   a_M: ratio of virtual cell height consisting of mobile CO2
    %   a_R: ratio of virtual cell height consisting of residual CO2 (NOT
    %        including residual part in mobile plume)
    %   sG: gas saturation

    z_RT = t + h_T; % bottom of residual part on TOP of col
    z_RB = t + h_B; % top of residual part on BOTTOM of col
    z_M = t + h; % mobile height
    H = B - T;
    
    a_M = max(min((z_M - T)./H, 1), 0); % scaling factor for plume saturation
    a_RB = max(min((B - z_RB)./H, 1), 0); % changed from B - h_B to B - z_RB
    a_RT = 1 - a_M - min(max((B - z_RT)./H, 0), 1); % scaling factor for residual saturation
    %a_RT = max((z_RT - max(T,z_M))./H, 0);
    
    a_R = max(a_RT + a_RB, 0).*(snr > 0); % max to avoid potential negative saturations
    sG = a_M.*(1 - swr) + a_R.*snr; % if below residual CO2 zone, sG = 0   
    
    if nargin > 9
        % --- OLD, DEPRECATED (?)---
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
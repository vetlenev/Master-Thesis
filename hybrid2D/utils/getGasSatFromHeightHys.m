function [a_M, a_R, sG] = getGasSatFromHeightHys(T, t, B, b, h_p, h_max, swr, snr, varargin)
    z_R = t + h_max;
    z_M = t + h_p;
    %z = (T + B)/2;
    H = B - T;
    
    a_M = max(min((z_M - T)./H, 1), 0); % scaling factor for plume saturation
    a_R = 1 - a_M - min(max((B - z_R)./H, 0), 1); % scaling factor for residual saturation
    
    sne = 0;
    sG = a_M.*(1 - swr) + a_R.*sne; % if below residual CO2 zone, sG = 0
    
    if nargin > 8
        % Specific reconstruction for NVE cells
        cNVE = varargin{1};            
        z_R = b(cNVE) - h_max(cNVE); % fill from bottom
        a_M(cNVE) = 0; % no mobile plume
        a_R(cNVE) = min(max((B(cNVE) - z_R)./H(cNVE), 0), 1);
        
        sG(cNVE) = a_R(cNVE).*snr;
    end
end
function [a_M, a_R, sG] = getGasSatFromHeight(T, t, B, b, h_p, h_max, swr, snr)
    z_R = t + h_max;
    z_M = t + h_p;
    %z = (T + B)/2;
    H = B - T;
    
    a_M = max(min((z_M - T)./H, 1), 0); % scaling factor for plume saturation
    a_R = 1 - a_M - min(max((B - z_R)./H, 0), 1); % scaling factor for residual saturation
    
    sG = a_M.*(1 - swr) + a_R.*snr; % if below residual CO2 zone, sG = 0
end
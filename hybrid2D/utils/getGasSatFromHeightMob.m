function [a_M, a_R, sG] = getGasSatFromHeightMob(T, t, B, b, h, h_T, h_B, swr, snr, varargin)
    z_RT = t + h_T; % bottom of residual part on TOP of col
    z_RB = t + h_B; % top of residual part on BOTTOM of col
    z_M = t + h; % mobile height
    H = B - T; % height of fine cell to reconstruct
    
    a_M = max(min((z_M - T)./H, 1), 0); % scaling factor for plume saturation
    a_RB = max(min((B - z_RB)./H, 1), 0); % changed from B - h_B to B - z_RB
    a_RT = 1 - a_M - min(max((B - z_RT)./H, 0), 1); % scaling factor for residual saturation
    
    a_R = max(a_RT + a_RB, 0).*(snr > 0); % max to avoid potential negative saturations
    sG = a_M.*(1 - swr) + a_R.*snr; % if below residual CO2 zone, sG = 0   
    
    if nargin > 9
        % Specific reconstruction for veHorizontal cells
        cHorz = varargin{1}; % veHorizontal transition cells satisfying t ~= T 
        cHorz_u = unique(cHorz);
        h_Bi = varargin{2}; % depth of (top of) bottom plume rising up from cell i
        H_i = varargin{3}; % depth of bottom of VE virtual cell i
        
        %a_RT = max(min((z_RT - T)./H, 1), 0);
        T = T(cHorz);
        B = B(cHorz);
        H = H(cHorz);
        
        a_M = a_M(cHorz_u);
        a_RB_top = accumarray(cHorz, min(max(h_Bi - T, 0), B-T)); % remove fraction at top of virtual cell not part of residual plume
        a_RB_top = a_RB_top(cHorz_u);
        a_RB_bottom = accumarray(cHorz, min(max(B - H_i, 0), B-T)); % remove corresponding part from bottom
        a_RB_bottom = a_RB_bottom(cHorz_u);
        a_RB = ((B - T) - a_RB_top - a_RB_bottom)./H;  
        
        a_RT = a_RT(cHorz_u);
        a_R = max(a_RT + a_RB, 0).*(snr > 0);
        
        sG(cHorz_u) = a_M.*(1-swr) + a_R.*snr;
        
        % Specific reconstruction for veBottom AND veHorizontal
        cBH = varargin{4};
        %cBH_h = cBH{1};  cBH_b = cBH{2};
        h_Bi = varargin{5};
        %h_Bh = h_Bi{1}; h_Bb = h_Bi{2};
        H_i = varargin{6};
        %H_ih = H_i{1}; H_ib = H_i{2}; % virtual VE cell height; VE cell height
        
        T = T(cBH);
        B = B(cBH);
        H = H(cBH);
        
        a_M = a_M(cBH);
        a_RB_top = min(max(h_Bi - T, 0), B-T); % remove fraction at top of virtual cell not part of residual plume
        %a_RB_top = a_RB_top(cHorz_u);
        a_RB_bottom = min(max(B - H_i, 0), B-T); % remove corresponding part from bottom
        %a_RB_bottom = a_RB_bottom(cHorz_u);
        a_RB = ((B - T) - a_RB_top - a_RB_bottom)./H;  
        
        a_RT = a_RT(cBH);
        a_R = max(a_RT + a_RB, 0).*(snr > 0);
        
        sG(cHorz_u) = a_M.*(1-swr) + a_R.*snr;
    end
end
function [a_Mob, a_R, sG] = getGasSatFromHeightFine(model, T, t, B, b, h, h_T, h_B, swr, snr, varargin)
    p = model.G.partition;
    
    z_RT = t + h_T; % bottom of residual part on TOP of col
    z_RB = t + h_B; % top of residual part on BOTTOM of col
    z_M = t + h; % mobile height
    H = B - T; % height of fine cell to reconstruct
    
    a_M = max(min((z_M - T)./H, 1), 0); % scaling factor for plume saturation
    a_RB = max(min((B - z_RB)./H, 1), 0); % changed from B - h_B to B - z_RB
    a_RT = 1 - a_M - min(max((B - z_RT)./H, 0), 1); % scaling factor for residual saturation
    
    a_R = max(a_RT + a_RB, 0).*(snr > 0); % max to avoid potential negative saturations    
    sG = a_M.*(1 - swr) + a_R.*snr; % residual part of mobile plume included in first term! 
    
    a_Mob = a_M;
    a_ResTop = a_RT;
    
    if nargin > 10        
        
        cHorz = varargin{1}; % veHorizontal transition cells satisfying t ~= T 
        cHorz_u = unique(cHorz);
        cBH = varargin{2};
        cBH_u = unique(cBH);
        
        hH_i = varargin{3}; % depth of (top of) bottom plume rising up from cell i
        hBH_i = varargin{4}; 
        
        H_i = varargin{5}; % depth of bottom of VE virtual cell i                                 
        BH_i = varargin{6};
               
        % --- Specific reconstruction for veHorizontal cells ---
        % --- Loop through each unique cHorz cell ---
        for c=1:numel(cHorz_u)
            ch = cHorz_u(c);
            pcH = p == ch; % fine cell index
            H_ik = (H_i(cHorz == ch))'; % transpose necessary fo doing copmutations for each component (virtual VE cell)
            hB_ik = (hH_i(cHorz == ch))';
            % -- Think this is ok --
            H_ik = min(t(pcH)) + H_ik;
            hB_ik = min(t(pcH)) + hB_ik;
            % ----------------------
            T_all = T(pcH);
            B_all = B(pcH);
                      
            partly_below = (T_all >= hB_ik & T_all < H_ik ...
                            & B_all > H_ik); % fine cell is partly below bottom residual plume
            partly_above = (T_all < hB_ik & ...
                            B_all <= H_ik & B_all > hB_ik); % fine cell partly above bottom residual plume
            outside = (T_all < hB_ik & B_all > H_ik); % fine cell covers entire residual plume and extends above and below it
            inside = (T_all >= hB_ik & T_all < H_ik ...
                        & B_all <= H_ik & B_all > hB_ik); % entire fine cell is inside residual plume
            
            a_RB = (H_ik - T_all).*partly_below + ...
                    (B_all - hB_ik).*partly_above + ...
                    (H_ik - hB_ik).*outside + ...
                    (B_all - T_all).*inside;  
                            
            a_RB = a_RB./(B_all-T_all); % scale by fine cell height to get in range [0,1]
            a_RB = sum(a_RB, 2); % sum rowwise to get one residual factor per fine cell
          
            a_M = a_Mob(pcH);
            a_RT = a_ResTop(pcH); 
            a_R_ = min(max(a_RT + a_RB, 0), 1).*(snr > 0);
        
            sG(pcH) = a_M.*(1-swr) + a_R_.*snr;
            % update residual factor, mobile factor same as before
            a_R(pcH) = a_R_;           
        end       
        
        % -- Specific reconstruction for veBottom AND veHorizontal ---        
        % --- Loop through each unique cBottomHorz ---
        for c=1:numel(cBH_u)
            cbh = cBH_u(c);
            pcH = p == cbh; % fine cell index
            H_ik = (BH_i(cBH == cbh))';
            hB_ik = (hBH_i(cBH == cbh))';
            % -- Test --
            H_ik = min(t(pcH)) + H_ik;
            hB_ik = min(t(pcH)) + hB_ik;
            % ----------
            T_all = T(pcH);
            B_all = B(pcH);
            
            partly_below = (T_all >= hB_ik & T_all < H_ik ...
                            & B_all > H_ik); % fine cell is partly below bottom residual plume
            partly_above = (T_all < hB_ik & ...
                            B_all <= H_ik & B_all > hB_ik); % fine cell partly above bottom residual plume
            outside = (T_all < hB_ik & B_all > H_ik); % fine cell covers entire residual plume and extends above and below it
            inside = (T_all >= hB_ik & T_all < H_ik ...
                        & B_all <= H_ik & B_all > hB_ik); % entire fine cell is inside residual plume
            
            a_RB = (H_ik - T_all).*partly_below + ...
                    (B_all - hB_ik).*partly_above + ...
                    (H_ik - hB_ik).*outside + ...
                    (B_all - T_all).*inside;  
                
            a_RB = a_RB./(B_all-T_all);
            a_RB = min(sum(a_RB, 2), 1); % sum rowwise to get one residual factor per fine cell, and cap at 1 if fine cell covers multiple residual regions
                      
            a_M = a_Mob(pcH);
            a_RT = a_ResTop(pcH);
            a_R_ = min(max(a_RT + a_RB, 0), 1-a_M).*(snr > 0);
            
            sG(pcH) = a_M.*(1-swr) + a_R_.*snr;
            a_R(pcH) = a_R_; % update residual content for relaxed VE cells
        end       
    end
end
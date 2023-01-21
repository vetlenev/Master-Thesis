function [a_Mob, a_R, sG] = getGasSatFromHeightFine_test(cells, T, t, B, b, h, h_T, h_B, swr, snr, varargin)
    %p = model.G.partition;   
    z_RT = t + h_T; % bottom of residual part on TOP of col
    z_RB = t + h_B; % top of residual part on BOTTOM of col
    z_M = t + h; % mobile height
    H = B - T; % height of fine cell to reconstruct
    
    a_M = max(min((z_M - T)./H, 1), 0); % scaling factor for plume saturation
    a_RB = max(min((B - z_RB)./H, 1), 0); % changed from B - h_B to B - z_RB
    a_RT = 1 - a_M - min(max((B - z_RT)./H, 0), 1); % scaling factor for residual saturation
    
    a_R = max(a_RT + a_RB, 0).*(snr > 0); % max to avoid potential negative saturations
    sG = a_M.*(1 - swr) + a_R.*snr; % if below residual CO2 zone, sG = 0   
    
    a_Mob = a_M;
    a_ResTop = a_RT;
    
    if nargin > 10        
        
        cHorz = varargin{1}; % veHorizontal transition cells satisfying t ~= T 
        cHorz_u = unique(cHorz);
        cBH = varargin{2};
        cBH_u = unique(cBH);
        
        hH_i = value(varargin{3}); % depth of (top of) bottom plume rising up from cell i
        hBH_i = value(varargin{4}); % cast from ADI to double to allow transposing
        
        H_i = varargin{5}; % depth of bottom of VE virtual cell i                                 
        BH_i = varargin{6};
               
        % --- Specific reconstruction for veHorizontal cells ---
        % --- Loop through each unique cHorz cell ---
        cHorz_nonempty = intersect(cHorz_u, cells); % only choose cHorz included in 'cells'
        
        for c=1:numel(cHorz_nonempty)
            ch = cHorz_nonempty(c);
            c_global = cells == ch;
            c_local = cHorz == ch;
%             if ~any(c_global) % none of input cells are relaxed VE cells => skip to next relaxed VE cell
%                continue; 
%             end
            H_ik = (H_i(c_local))'; % transpose necessary fo doing copmutations for each component (virtual VE cell)
            hB_ik = (hH_i(c_local))';
            % -- Think this is ok --
            H_ik = min(t(c_global)) + H_ik;
            hB_ik = min(t(c_global)) + hB_ik;
            % ----------------------
            T_all = T(c_global);
            B_all = B(c_global);
                      
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
          
            a_M = a_Mob(c_global);
            a_RT = a_ResTop(c_global); 
            a_R_ = min(max(a_RT + a_RB, 0), 1).*(snr > 0);
        
            sG(c_global) = a_M.*(1-swr) + a_R_.*snr;
            % update residual factor, mobile factor same as before
            a_R(c_global) = a_R_;           
        end       
        
        % -- Specific reconstruction for veBottom AND veHorizontal ---        
        % --- Loop through each unique cBottomHorz ---
        cBH_nonempty = intersect(cBH_u, cells); % only choose cBH included in 'cells'
        
        for c=1:numel(cBH_nonempty)
            cbh = cBH_nonempty(c);
            c_global = cells == cbh;
            c_local = cBH == cbh;
%             if ~any(c_global) % none of input cells are relaxed VE cells => skip to next relaxed VE cell
%                continue; 
%             end
            %pcH = p == cbh; % fine cell index
            H_ik = (BH_i(c_local))';
            hB_ik = (hBH_i(c_local))';
            % -- Test --
            H_ik = min(t(c_global)) + H_ik; % height from top of parent cell down to sealing layer
            hB_ik = min(t(c_global)) + hB_ik; % height from top of parent cell to top of plume
            % ----------
            T_all = T(c_global);
            B_all = B(c_global);
            
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
                      
            a_M = a_Mob(c_global);
            a_RT = a_ResTop(c_global);
            a_R_ = min(max(a_RT + a_RB, 0), 1-a_M).*(snr > 0);
            
            sG(c_global) = a_M.*(1-swr) + a_R_.*snr;
            a_R(c_global) = a_R_; % update residual content for relaxed VE cells
        end       
    end
end
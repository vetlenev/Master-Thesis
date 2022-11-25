function [a_Mob, a_R, sG] = getGasSatFromHeightMob(model, T, t, B, b, h, h_T, h_B, swr, snr, varargin)
    p = model.G.partition;
    
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
            
            %a_RBi = zeros(nnz(pcH), 1);
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
                
            %no_intersect = (T_all >= H_ik | B_all <= hB_ik);
            %a_RB = a_RB.*(~no_intersect); % zero out fine cells not intersecting bottom residual plume
            a_RB = a_RB./(B_all-T_all); % scale by fine cell height to get in range [0,1]
            a_RB = sum(a_RB, 2); % sum rowwise to get one residual factor per fine cell
            
            %a_RB_top = sum(min(max(B_all-H_ik', 0), B_all-T_all), 2); % transpose of H_ik to get array of plume heights for each fine cell
            %a_RB_bottom = sum(min(max(hB_ik'-T_all, 0), B_all-T_all), 2); % transpose here as well
            %a_RB = (B_all-T_all - a_RB_top - a_RB_bottom)./H(pcH);
            a_M = a_Mob(pcH);
            a_RT = a_ResTop(pcH); 
            a_R_ = min(max(a_RT + a_RB, 0), 1).*(snr > 0);
        
            sG(pcH) = a_M.*(1-swr) + a_R_.*snr;
            % update residual factor, mobile factor same as before
            a_R(pcH) = a_R_;           
        end
        % -------------------------------------------
        
%         T_all = T(cHorz); T_u = T(cHorz_u);
%         B_all = B(cHorz); B_u = B(cHorz_u);       
%         H_all = H(cHorz); H_u = H(cHorz_u);
%         
%         a_M = a_Mob(cHorz_u);
%         a_RB_top = accumarray(cHorz, min(max(hB_i - T_all, 0), B_all - T_all)); % remove fraction at top of virtual cell not part of residual plume
%         a_RB_top = a_RB_top(cHorz_u);
%         a_RB_bottom = accumarray(cHorz, min(max(B_all - H_i, 0), B_all - T_all)); % remove corresponding part from bottom
%         a_RB_bottom = a_RB_bottom(cHorz_u);        
%         a_RB = ((B_u - T_u) - a_RB_top - a_RB_bottom)./H_u;  
%         
%         a_RT = a_ResTop(cHorz_u);
%         a_R = max(a_RT + a_RB, 0).*(snr > 0);
%         
%         sG(cHorz_u) = a_M.*(1-swr) + a_R.*snr;
        % -----------------
               
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
                
            %no_intersect = (T_all >= H_ik | B_all <= hB_ik);
            %a_RB = a_RB.*(~no_intersect); % zero out fine cells not intersecting bottom residual plume
            a_RB = a_RB./(B_all-T_all);
            a_RB = min(sum(a_RB, 2), 1); % sum rowwise to get one residual factor per fine cell, and cap at 1 if fine cell covers multiple residual regions
                      
            a_M = a_Mob(pcH);
            a_RT = a_ResTop(pcH);
            a_R_ = min(max(a_RT + a_RB, 0), 1-a_M).*(snr > 0);
            
            %a_R_(a_M == 1) = 0; % force zero residual part in fully mobile zone
        
            sG(pcH) = a_M.*(1-swr) + a_R_.*snr;
            a_R(pcH) = a_R_;
        end
        % --------------------------------------------
       
%         T = T(cBH);
%         B = B(cBH);
%         H = H(cBH);
%         
%         a_M = a_Mob(cBH);
%         a_RB_top = min(max(hBH_i - T, 0), B-T); % remove fraction at top of virtual cell not part of residual plume
%         %a_RB_top = a_RB_top(cHorz_u);
%         a_RB_bottom = min(max(B - BH_i, 0), B-T); % remove corresponding part from bottom
%         %a_RB_bottom = a_RB_bottom(cHorz_u);
%         a_RB = ((B - T) - a_RB_top - a_RB_bottom)./H;  
%         
%         a_RT = a_ResTop(cBH);
%         a_R = max(a_RT + a_RB, 0).*(snr > 0);
%         
%         sG(cBH) = a_M.*(1-swr) + a_R.*snr;
        % --------------------
    end
end
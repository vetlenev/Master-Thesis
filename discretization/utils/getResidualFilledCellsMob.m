function [c_NVE, c_NVE_Mob, c_VE] = getResidualFilledCellsMob(model, pv, bG, Sn, Sn_max, vG_sum)
    % Get VE columns that are filled with CO2 from bottom.
    % Used for correct computation of fluxes in the residual framework.
    %
    % Inputs:
    %   T/B: top/bottom of virtual cells
    %   t/b: top/bottom of column
    %   h: global height of plume in column
    %   h_max: global max depth CO2 has reached
    %   swr/snr: residual water/co2 saturation
    %   veBottomNonzeroSn: boolean array indicating cells in vertical
    %                      ve transition zone, where ve is upper part and ve/fine cells lower part.
    %   Sn: Co2 saturation from COARSE grid   
    %   vG_sum [m^3]: summed CO2 flux (abs value) over time from bottom interface 
    %                up to current time step  
    
    CG = model.G;
    op = model.operators;    
    isFine = CG.cells.discretization == 1;      
    
    p = CG.partition;
    %Sn_p = Sn(p);
    %snMax_p = snMax(p);
    
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
    %[a_M, a_R, sG] = getGasSatFromHeight(TT, tt, BB, bb, h, h_max, swr, snr);
    % ------------------- 
      
    % Modify for ve bottom transition
    veTransition = op.connections.veToFineVertical | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;
    veInternalConn = op.connections.veInternalConn; 
    
    c_NVE = {};
    c_NVE_Mob = {};
    c_VE = {};
    
    n = op.N; 
    %n = CG.faces.neighbors;
    cn = op.N(veTransition, :);
    vG_bottom = abs(vG_sum(ismember(n, cn, 'rows'))); % take abs val since fluxes are negative from bottom and up
    
    % NEW
    for idx=1:2 % CONSIDER REMOVING LOOP: veBottomFilled and veTopEmpty will never both occur as it stand now!
        c = cn(:, idx);        
        isVE_c = ~isFine(c);
           
        if any(isVE_c) % reconstruct specifically for ve col transitioning to a lower layer
            t = CG.cells.topDepth(c);
            T = op.connections.faceTopDepth(veTransition, idx);
            b = CG.cells.bottomDepth(c);
            B = op.connections.faceBottomDepth(veTransition, idx);
            H = CG.cells.height(c);
            
            % Compare pore volume and darcy flux of co2
            S = value(Sn(c)); 
            S_max = value(Sn_max(c));            
            S_B = 1./pv(c) .*vG_bottom; % CO2 sat originating from bottom flux
            S_max_T = S_max - S_B; % max saturation reached for top part (excluding part originating from bottom flux)
            
            h_T = H.*S_max_T./(1-swr);
            h_B = max(h_T, H.*(1-S_B./snr));
            
            Snr_tot = snr./H.*(h_T + (H - h_B)); % total residual saturation in column (NB: includes residual part in mobile zone!)
            S_mob = S - Snr_tot;
            
            h = H.*(S_mob./(1-swr-snr)); % subtract snr since this was included in Snr_tot
          
            buff = 0.05;
            % T =~ t to avoid potentially selecting fine columns where B ==
            % b and T == t, but T =~ t doesn't matter for VE columns
            veBottom = (B == b & T ~= t) & (vG_bottom > 0 & S < snr & vG_bottom > (1-buff)*co2_vol); % & vG_veTrans < (1+buff)*co2_vol); % co2 ONLY filled from bottom, and has not yet reached top => not in VE!                      
            % Check if mobile top plume has reached residual bottom part
            % => VE holds!            
            veValid = (B == b & T ~= t) & (vG_bottom > 0 & h_B + h_T >= (1-buff).*H); % residual plume has reached top -> from prf to frf -> VE assumption holds
                                             
            %c_veNot = intersect(c(veBottomFilled), c(veTopEmpty)); % internal ve columns where ve assumption not satisfied
            c_nve = c(veBottom);
            c_nve_mob = c(veBottomTop); % NVE cells with finite mobile plume on top
            c_ve = c(veValid);
            
            c_vic = op.N(veInternalConn, :);            
            
            c_nve = c_vic(ismember(c_vic, c_nve));           
            if ~isempty(c_nve)
                c_NVE = cat(2, c_NVE, c_nve);
            end
            
            c_nve_mob = c_vic(ismember(c_vic, c_nve_mob));           
            if ~isempty(c_nve_mob)
                c_NVE_Mob = cat(2, c_NVE_Mob, c_nve_mob);
            end
            
            c_ve = c_vic(ismember(c_vic, c_ve));
            if ~isempty(c_ve)
                c_VE = cat(2, c_VE, c_ve);                            
            end
        end
    end
    
    % OLD
%     for idx=1:2 % CONSIDER REMOVING LOOP: veBottomFilled and veTopEmpty will never both occur as it stand now!
%         c = cn(:, idx);        
%         isVE_c = ~isFine(c);
%            
%         if any(isVE_c) % reconstruct specifically for ve col transitioning to a lower layer
%             t = CG.cells.topDepth(c);
%             T = op.connections.faceTopDepth(veTransition, idx);
%             b = CG.cells.bottomDepth(c);
%             B = op.connections.faceBottomDepth(veTransition, idx);
%             H = CG.cells.height(c);
%             
%             Compare pore volume and darcy flux of co2
%             Sn_c = value(Sn(c));
%             pv = poreVolume(CG, model.rock);
%             pv = pv(c);
%             pv_bG = pv.*bG(c);
%             co2_vol = pv_bG.*Sn_c;
%             
%             h_B = min(vG_bottom./(pv_bG.*snr).*H, H); % bottom residual part
%             h_T = max((co2_vol - vG_bottom)./(pv_bG.*(1-swr)).*H, 0); % top mobile part
%             
%             buff = 0.05;
%             veBottom = (B == b & T ~= t) & (vG_bottom > 0 & Sn_c < snr & vG_bottom > (1-buff)*co2_vol); % & vG_veTrans < (1+buff)*co2_vol); % co2 ONLY filled from bottom, and has not yet reached top => not in VE!          
%             veBottomTop = (B == b & T ~= t) & (vG_bottom > 0 & vG_bottom <= (1-buff)*co2_vol ...
%                                                               & h_B + h_T < (1-buff).*H);
%             
%             Check if mobile top plume has reached residual bottom part
%             => VE holds!            
%             veValid = (B == b & T ~= t) & (vG_bottom > 0 & h_B + h_T >= (1-buff).*H); % residual plume has reached top -> from prf to frf -> VE assumption holds
%                                              
%             c_veNot = intersect(c(veBottomFilled), c(veTopEmpty)); % internal ve columns where ve assumption not satisfied
%             c_nve = c(veBottom);
%             c_nve_mob = c(veBottomTop); % NVE cells with finite mobile plume on top
%             c_ve = c(veValid);
%             
%             c_vic = op.N(veInternalConn, :);            
%             
%             c_nve = c_vic(ismember(c_vic, c_nve));           
%             if ~isempty(c_nve)
%                 c_NVE = cat(2, c_NVE, c_nve);
%             end
%             
%             c_nve_mob = c_vic(ismember(c_vic, c_nve_mob));           
%             if ~isempty(c_nve_mob)
%                 c_NVE_Mob = cat(2, c_NVE_Mob, c_nve_mob);
%             end
%             
%             c_ve = c_vic(ismember(c_vic, c_ve));
%             if ~isempty(c_ve)
%                 c_VE = cat(2, c_VE, c_ve);                            
%             end
%         end
%     end
    
    c_NVE = unique(cell2mat(c_NVE));
    c_NVE_Mob = unique(cell2mat(c_NVE_Mob));
    c_VE = unique(cell2mat(c_VE));             
       
end

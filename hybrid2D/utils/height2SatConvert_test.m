function [sG] = height2SatConvert_test(model, h_func, h_max_func, h_T, h_B, Sn, snMax, vG, ii)
    % Transformation from from plume and residual depths to coarse
    % saturations. Used for convertion from coarse to fine states.
    % Inputs:
    %   T/B: top/bottom of virtual cells
    %   t/b: top/bottom of column
    %   h: global height of plume in column
    %   h_max: global max depth CO2 has reached
    %   swr/snr: residual water/co2 saturation
    %   veBottomNonzeroSn: boolean array indicating cells in vertical
    %                      ve transition zone, where ve is upper part and ve/fine cells lower part.
    %   Sn: Co2 saturation from COARSE grid
    CG = model.G;
    op = model.operators;
    isFine = CG.cells.discretization == 1;      
    
    p = CG.partition;
    Sn_p = Sn(p);
    snMax_p = snMax(p);  
    
    tt = CG.cells.topDepth(p);
    TT = CG.parent.cells.topDepth;

    bb = CG.cells.bottomDepth(p);
    BB = CG.parent.cells.bottomDepth;
    H = CG.cells.height(p);
    
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
    
    % --- Standard ---    
    %h = h_func(Sn_p, snMax_p, H);
    %h_max = h_max_func(snMax_p, H);
    h = h_func(p);
    h_T = h_T(p);
    h_B = h_B(p);
    
    % Works for all cells -> h_B = H for no bottom flux
    [a_M, a_R, sG] = getGasSatFromHeightMob2(TT, tt, BB, bb, h, h_T, h_B, swr, snr);
    % -------------------
      
    %[c_VE_Not, c_VE, c_VE_Horz] = getResidualFilledCells_test(model, Sn, vG, Sn0, sgNVE0); % coarse saturation Sn, not fine saturation sG
    c_VE_Not = []; c_VE = []; c_VE_Horz = [];
    
    if ii == 220
       test = 0; 
    end
       
    if ~isempty(c_VE_Not)
        % --- 1. Partly residual filled (prf) ---               
        partly_res_filled = c_VE_Not; % fill from BOTTOM
        
        p_prf = ismember(p, partly_res_filled); % fine-scale cells that are partly filled     
               
        t_p = tt(p_prf); % fine-scale global ve heights        
        T_p = TT(p_prf); % fine-scale local cell heights        
        b_p = bb(p_prf);       
        B_p = BB(p_prf);
        H = b_p - t_p; % height of ve column partly filled      
        
        h_prf = H .* Sn_p(p_prf)./snr; % ratio of column height filled with residual co2 from below        
        % NB: Special treatment of sG!   
        sG_prf = fillSatFromBelow(t_p, T_p, b_p, B_p, h_prf, snr); % sGmax not needed as input; it is forced to 1
               
        sG(p_prf) = sG_prf; % update fine-scale saturation for partly filled cells      
        % ---------------------------------
    end
    
    if ~isempty(c_VE) % only update reconstructed saturation if co2 originates from bottom of any ve column
        % --- 2. Fully residual filled (frf) ---                
        fully_res_filled = c_VE; % fill from TOP       
        
        p_frf = ismember(p, fully_res_filled);
              
        % T, t, B and b defined for coarse cells
        t_f = tt(p_frf);           
        T_f = TT(p_frf);
        b_f = bb(p_frf);
        B_f = BB(p_frf);
        H = b_f - t_f;                
        
        h_frf = h_func(Sn_p(p_frf), 1-swr, H);
        h_max_frf = h_max_func(1-swr, H); % plume originates from bottom -> enforce snGmax = 1-swr 

        [~, ~, sG_frf] = getGasSatFromHeight(T_f, t_f, B_f, b_f, h_frf, ...
                                            h_max_frf, swr, snr);
     
        sG(p_frf) = sG_frf; % update fine-scale saturation for fully filled cells     
        % -------------------------------- 
    end
    
%     if ~isempty(c_VE_Mob)
%        % --- 3. NVE cells with finite mobile plume at top --- 
%     end
    
    if ~isempty(c_VE_Horz)
        % --- 4. NVE cells with horizontal fluxes from neighbors ---
        p_nve_horz = ismember(p, c_VE_Horz);
        
        t_h = tt(p_nve_horz);           
        T_h = TT(p_nve_horz);
        b_h = bb(p_nve_horz);
        B_h = BB(p_nve_horz);
        H = b_h - t_h; 
        
        Sn0_h = Sn0_p(p_nve_horz);
        sGnve0_h = sgNVE0_p(p_nve_horz);
        
        h = (Sn0_h - sGnve0_h).*H./(1-swr);
        h_max = sGnve0_h.*H./snr + h;
        
        sG_nve_horz = fillNVEHorizontal(t_h, T_h, b_h, B_h, h, h_max, swr, snr);
        sG(p_nve_horz) = sG_nve_horz;
    end
    
end

function [sG] = fillSatFromBelow(t, T, b, B, h, snr)
    %Sn = Sn(vePartlyFilled);
    z_R = b - h; % top of immobile plume
    ratio_inside_cell = ((B - z_R)./(B - T));
    sG = snr.*(T >= z_R) + snr.*ratio_inside_cell.*(T < z_R & z_R < B);
end

function [sG] = fillNVEHorizontal(t, T, b, B, h, h_max, swr, snr)
    sG = zeros(numel(h), 1);
    z_P = t + h; % bottom of plume
    z_R = b - (h_max - h); % top of residual part
    
    sG(T > z_R) = snr; % entire cell in residual part
    sG(B < z_P) = 1-swr; % entire cell in mobile plume
        
    ratio_R = 1 - (z_R - T)./(B - T);
    ratio_P = (z_P - T)./(B - T);
    
    sG(T <= z_R & B >= z_R) = ratio_R.*snr;
    sG(T <= z_P & B >= z_R) = ratio_P.*(1-swr);
end

function [sG, h, h_max, c_VE_Not] = height2SatConvert(model, h_func, h_max_func, Sn, snMax, vG, ii)
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
    % Returns:
    %   sG: fine-scale gas saturation reconstructed from height
    %   h: global height of plume for fine cells
    %   h_max: global max height of plume through history
    %   c_VE_Not: cells satisfying NVE
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
    
    [c_VE_Not, c_VE] = getResidualFilledCells(model, Sn, vG); % coarse saturation Sn, not fine saturation sG
    cNVE = ismember(p, c_VE_Not); % get parent cells
    
    % --- 3. Standard ---    
    h = h_func(Sn_p, snMax_p, H);
    h_max = h_max_func(snMax_p, H);    
    
    [a_M, a_R, sG] = getGasSatFromHeight(TT, tt, BB, bb, h, h_max, swr, snr, cNVE);
    % -------------------
          
    if ii == 140
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
        % ---
        h(p_prf) = 0;
        h_max(p_prf) = h_prf;
        % ---
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
        % ---
        h(p_frf) = h_frf;
        h_max(p_frf) = h_max_frf;
        % ---

        [~, ~, sG_frf] = getGasSatFromHeight(T_f, t_f, B_f, b_f, h_frf, ...
                                            h_max_frf, swr, snr);
     
        sG(p_frf) = sG_frf; % update fine-scale saturation for fully filled cells     
        % -------------------------------- 
    end       
    
end

function [sG] = fillSatFromBelow(t, T, b, B, h, snr)
    %Sn = Sn(vePartlyFilled);
    z_R = b - h; % top of immobile plume
    ratio_inside_cell = ((B - z_R)./(B - T));
    sG = snr.*(T >= z_R) + snr.*ratio_inside_cell.*(T < z_R & z_R < B);
end

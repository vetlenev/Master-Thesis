function [sG, snMaxVE_bottom, vG_any] = height2SatConvert(model, h_func, h_max_func, swr, snr, Sn, ...
                                            snMax, snMaxVE_bottom, vG, vG_any, ii)
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
    
    % --- 3. Standard ---    
    h = h_func(Sn_p, snMax_p, H);
    h_max = h_max_func(snMax_p, H);
    
    [a_M, a_R, sG] = getGasSatFromHeight(TT, tt, BB, bb, h, h_max, swr, snr);
    % -------------------
      
    % Modify for ve bottom transition
    veTransition = op.connections.veToFineConn | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;
    veInternalConn = op.connections.veInternalConn; 
      
    %n = CG.faces.neighbors;
    n = op.N;
    cn = op.N(veTransition, :);
    vG_veTrans = vG(ismember(n, cn, 'rows'));
    
    %vG_any = max(vG_veTrans > 0, vG_any); % if any vG goes from zero to non-zero, vG_any updated to 1
    
    for idx=1:2 % loop over neighbor set
        c = cn(:, idx);         
        isVE_c = ~isFine(c);
        
        if any(isVE_c) % reconstruct specifically for ve col transitioning to a lower layer
            t = CG.cells.topDepth(c);
            T = op.connections.faceTopDepth(veTransition, idx);
            b = CG.cells.bottomDepth(c);
            B = op.connections.faceBottomDepth(veTransition, idx);                                                             
                                  
            veBottomConnIdx = (B == b & T ~= t) & (vG_veTrans > 0 & Sn_c < snr); % transition from UPPER ve column to LOWER fine cells, not from LOWER ve col to UPPER fine cells, and plume comes from lower layer.
            %veBottomConnIdx = (B == b & T ~= t) & (abs(vG_veTrans) > 0 | vG_any); % nonzero flux of co2 from bottom
            veValid = (B == b & T ~= t) & (vG_veTrans > 0 & Sn_c >= snr);
            
            c_veNot = c(veBottomConnIdx);
            c_ve = v(veValid);
   
            c_vic = op.N(veInternalConn, :);
            c_veNot_bool = ismember(c_vic, c_veNot);
            c_ve_bool = ismember(c_vic, c_ve);
            c_VE_Not = c_vic(c_veNot_bool);
            c_VE = c_vic(c_ve_bool);
                 
        end
    end
    
    c_VE_Not = unique(c_VE_Not);
    
    if ii == 140
       test = 0; 
    end
       
    if ~isempty(c_VE_Not) % only update reconstructed saturation if co2 originates from bottom of any ve column
        % --- 2. Fully residual filled (frf) ---
        c_vic = unique(c_vic); % all unique internal ve cells
        c_VE_Not = unique(c_VE_Not); % all unique internal ve cells whose bottom is filled with any CO2
        Sn_vic = Sn(c_vic);
        Sn_vic(~ismember(c_vic, c_VE_Not)) = 0; % don't include internal VE cells whose plume doesn't originate from below
        
        snMaxVE_bottom = max(Sn_vic, snMaxVE_bottom); % max residual part obtained so far
        snMaxVE_bottom_cve = snMaxVE_bottom(c_VE_Not);
       
        frf = snMaxVE_bottom_cve >= snr;
        fully_res_filled = c_VE_Not(frf); % fill from TOP       
        
        %fully_res_filled = c_VE(Sn(c_VE) >= snr);
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

        % NB: sG defined for parent grid, transform: coarse -> fine          
        sG(p_frf) = sG_frf;      
        % --------------------------------          
    
        % --- 1. Partly residual filled (prf) ---       
        prf = snMaxVE_bottom_cve < snr;
        partly_res_filled = c_VE_Not(prf); % fill from BOTTOM       
        
        p_prf = ismember(p, partly_res_filled); % fine-scale cells that are partly filled     
               
        t_p = tt(p_prf); % fine-scale global ve heights        
        T_p = TT(p_prf); % fine-scale local cell heights        
        b_p = bb(p_prf);       
        B_p = BB(p_prf);
        H = b_p - t_p; % height of ve column partly filled      
        
        h_prf = H .* Sn_p(p_prf)./snr; % ratio of column height filled with residual co2 from below        
        % NB: Special treatment of sG!   
        sG_prf = fillSatFromBelow(t_p, T_p, b_p, B_p, h_prf, snr); % sGmax not needed as input; it is forced to 1
               
        sG(p_prf) = sG_prf;       
        % ---------------------------------
    end
    
end

function [sG] = fillSatFromBelow(t, T, b, B, h, snr)
    %Sn = Sn(vePartlyFilled);
    z_R = b - h; % top of immobile plume
    ratio_inside_cell = ((B - z_R)./(B - T));
    sG = snr.*(T >= z_R) + snr.*ratio_inside_cell.*(T < z_R & z_R < B);
end

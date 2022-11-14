function [vW, vG, mobW, mobG, upcw, upcg, h_global] = computeHybridFluxesVEpc(model, pW, sG, muW, muG, rhoW, rhoG, trans, sgMax, cellsVENot, cellsVE, cellsVEHorz, varargin)
% Internal function - computes interior fluxes for the hybrid VE models

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
    opt = struct('ve_model', 'sharp_interface'); % use sharp interface by default
    opt = merge_options(opt, varargin{:});
    
    op = model.operators;
    G = model.G;
    f = model.fluid;
    g = norm(model.gravity);   
    
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
    isFine = G.cells.discretization == 1;  
    all_coarse_cells = (1:G.cells.num)';
    
    sgMax_dummy = sgMax; % used to assign correct h_max to VE cells filled from bottom
    
    if ~isempty(cellsVE) % once residual plume reaches top, set sgMax to maximum and use standard formulas       
        sgMax_dummy(cellsVE) = 1-swr; % would otherwise be discontinuous transition from H to snr*H
    end 
            
    % --- upscaledSat2height ---
    if strcmp(opt.ve_model, 'sharp_interface')
        % Sharp interface: Plume height simple expression of sG
        H = G.cells.height;        
        h_max_global = H.*sgMax_dummy./(1-swr); 
        % height of mobile plume has special formula -> takes into account
        % CO2 sat of VE column and max CO2 sat (i.e. need to SUBTRACT
        % residual part from current net CO2 sat)
        h_global = H.*((sG.*(1-swr) - sgMax_dummy.*snr) ./ ((1-swr).*(1-swr-snr)));    
    elseif strcmp(opt.ve_model, 'general')
        % General model: Plume height evaluated from upscaled pc:
        % Pc = drho * g * h
        assert(~isempty(f.pcWG));
        pc    = f.pcWG(sG, pW, 'sGmax', sgMax_dummy);
        pcMax = f.pcWG(sgMax_dummy, pW, 'sGmax', sgMax_dummy);
        g_drho  = g * (rhoW(pW) - rhoG(pW));
        h_global     = pc ./ g_drho;
        h_max_global = pcMax ./ g_drho;
    else
        error('VE model %s not supported', opt.ve_model);
    end
    
    if ~isa(h_global, 'ADI')
       stop = 0; 
    end
    
    if ~isempty(cellsVENot) % only change height for transition zones
       % NB: Coarse saturation remains unchanged, but the DISTRIBUTION
       % (represented by h and h_max) is changed => modifies the fluxes
       h_global(cellsVENot) = 0; % no mobile part of plume => no CO2 mobility!
       %h_max_global(cellsVENot) = value(sG(cellsVENot)).*H(cellsVENot)./snr; % distribute "removed" mobile part to residual content to preserve mass 
       h_max_global(cellsVENot) = sgMax_dummy(cellsVENot).*H(cellsVENot)./snr;
    end     
    
    if ~isempty(cellsVE) && model.fluid.hys
       % h_max_global remains unchanged -> filled to bottom by default since sgMax = 1-swr
       C = 1/snr - 1/(1-swr);
       sne = sgMax(cellsVE) ./ (C*sgMax(cellsVE) + 1); % NB: use sgMax not sgMax_dummy since we want the TRUE max sat, not 1-swr
       h_global(cellsVE) = max(H(cellsVE).*((sG(cellsVE) - sne) ./ (1-swr-snr)), 0); 
    end
    % --------------------------               
    
    p_entry = f.pe_rest;
    
    %[pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h, H, rhoW, rhoG, muW, muG, isFine, isFine);
    [pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, pW, sG, h_global, h_max_global, H, rhoW, rhoG, muW, muG, ...
                                                isFine, all_coarse_cells, h_global, h_max_global, {}, {}, ...
                                                cellsVENot, 've_model', opt.ve_model, 'p_entry', p_entry);
    if isa(model, 'ThreePhaseCompositionalModel')
        sW = 1-sG;
        rhoWf = op.faceAvg(rhoW.*sW)./max(op.faceAvg(sW), 1e-8);
        rhoGf = op.faceAvg(rhoG.*sG)./max(op.faceAvg(sG), 1e-8);
    else
        rhoWf = op.faceAvg(rhoW);
        rhoGf = op.faceAvg(rhoG);
    end

    z = (G.cells.bottomDepth + G.cells.topDepth)/2;
    gdz =  g.*op.Grad(z);
    
    vIc = op.connections.veInternalConn;
    n1 = op.N(vIc, 1);
    n2 = op.N(vIc, 2);
     
    [gdz_g, gdz_w] = deal(gdz);
    % Compute upscaled height difference for internal ve columns, at reference pressure heights for CO2 and water:
    gdz_g(vIc) = g*(G.cells.topDepth(n2) - G.cells.topDepth(n1));
    % if water pressure computed at bottom of VE cells
    gdz_w(vIc) = g*(G.cells.bottomDepth(n2) - G.cells.bottomDepth(n1));        
    % if water pressure computed at top of VE cells
    %gdz_w(vIc) = g*(G.cells.topDepth(n2) - G.cells.topDepth(n1));
    
    dpG   = op.Grad(pG) - rhoGf .* gdz_g;
    upcg  = (value(dpG)<=0);
    vG = -op.faceUpstr(upcg, mobG) .* trans .* dpG;

    dpW   = op.Grad(pW) - rhoWf .* gdz_w;
    upcw  = (value(dpW)<=0);
    vW = -op.faceUpstr(upcw, mobW) .* trans .* dpW;

    % Treat transition between different discretizations
    transition_ve = model.operators.connections.veTransitionHorizontalConn & ~model.operators.connections.veToFineConn;
     
    if any(transition_ve)
        [vW(transition_ve), vG(transition_ve),...
         upcw(transition_ve), upcg(transition_ve)] = computeTransitionFluxVE(model, pW, h_global, h_max_global, rhoW, rhoG, muW, muG, ...
                                                                              transition_ve, true, cellsVENot, 'p_entry', p_entry);
    end
    % Treat transition between ve zones in vertical direction (for diffuse
    % leakage)
    transition_vertical = model.operators.connections.veToFineConn | ...
                          model.operators.connections.veTransitionVerticalConn & op.T > 0;
  
    if any(transition_vertical)        
        [vW(transition_vertical), vG(transition_vertical),...
         upcw(transition_vertical), upcg(transition_vertical)] = computeTransitionFluxVE(model, pW, h_global, h_max_global, rhoW, rhoG, muW, muG, ...
                                                                                            transition_vertical, false, cellsVENot, 'p_entry', p_entry);
    end

end

function [pW, pW_f, sG, h, h_max, H, rhow, rhog, muw, mug, isFine] = getTransitionValuesVE_coarse(model, pW, h_global, h_max_global, index, subs, rhoW, rhoG, muW, muG, cellsNVE, varargin)    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0);
    opt = merge_options(opt, varargin{:});
    
    c = model.operators.N(subs, index);
    cNVE = ismember(c, cellsNVE);
    
    isFine = model.G.cells.discretization(c) == 1;

    t = model.G.cells.topDepth(c);
    T = model.operators.connections.faceTopDepth(subs, index);
    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);

    %sG = getGasSaturationFromHeight(T, t, B, b, h_global(c));   
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
    [a_M, a_R, sG] = getGasSatFromHeight(T, t, B, b, h_global(c), h_max_global(c), swr, snr, cNVE); 
    
    H = model.operators.connections.faceHeight(subs, index);
    h = a_M.*H;  
    h_max = (a_R + a_M).*H;
    
    g = norm(model.gravity);
    rhow = rhoW(c);
    rhog = rhoG(c);
    
    pW_f = pW(c);
    muw = muW(c);
    mug = muG(c);      
    
    %pW = getWaterPressureFromHeight(B, t, B, b, h_global(c), h_max_global(c), pW(c), swr, snr, g, rhow, rhog, 'cNVE', cNVE, 'p_entry', opt.p_entry);   
    %pW = getPwFromHeight(B, t, b, h_global(c), h_max_global(c), pW(c), g, rhow, rhog, swr, snr, cellsNVE, c);
    pW = getPwFromHeight(B, t, b, h_global(c), h_max_global(c), pW(c), g, rhow, rhog, swr, snr, ...
                            'cNVE', cNVE, 'p_entry', opt.p_entry, 'res_type', []);
end

function [pW, pW_f, sG, h, h_max, H, rhow, rhog, muw, mug, isFine] = getTransitionValuesVE_fine(model, pW, h_global, h_max_global, index, subs, rhoW, rhoG, muW, muG, cellsNVE, varargin)    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0);
    opt = merge_options(opt, varargin{:});
    
    c = model.operators.N(subs, index);
    cNVE = ismember(c, cellsNVE);    
    t = model.G.cells.topDepth(c);
    
    isFine = model.G.cells.discretization(c) == 1;

    T = model.operators.connections.faceTopDepth(subs, index);

    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);
    H = model.operators.connections.faceHeight(subs, index);
    
    %sG = getGasSaturationFromHeight(T, t, B, b, h_global(c));
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);    
       
    [a_M, a_R, sG] = getGasSatFromHeight(T, t, B, b, h_global(c), h_max_global(c), swr, snr, cNVE);   
    
    h = a_M.*H;  
    h_max = (a_R + a_M).*H;
    
    g = norm(model.gravity);
    rhow = rhoW(c);
    rhog = rhoG(c);
    C = (T + B)/2;
    pW_f = pW(c);
    muw = muW(c);
    mug = muG(c);    
   
    %pW_c = getWaterPressureFromHeight(C, t, B, b, h_global(c), h_max_global(c), pW_f, swr, snr, g, rhow, rhog, 'cNVE', cNVE, 'p_entry', opt.p_entry);   
    %pW_c = getPwFromHeight(C, t, b, h_global(c), h_max_global(c), pW_f, g, rhow, rhog, swr, snr, cellsNVE, c);
    pW_c = getPwFromHeight(C, t, b, h_global(c), h_max_global(c), pW_f, g, rhow, rhog, swr, snr, ...
                            'cNVE', cNVE, 'p_entry', opt.p_entry, 'res_type', []);
    
    pW = isFine.*pW_f + ~isFine.*pW_c;
end


function [vW, vG, upcw, upcg] = computeTransitionFluxVE(model, pW, h, h_max, rhoW, rhoG, muW, muG, vtc, treatAsCoarse, cellsVENot, varargin)
    % INPUTS:
    %   h: global depth of mobile plume
    %   h_max: global height of residual part
    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0);
    opt = merge_options(opt, varargin{:});

    op = model.operators;
    G = model.G;
    trans = model.operators.T(vtc);
    
    % Consistent with rest of MRST...
    grad = @(l, r) r - l;
    upstr = @(flag, l, r) flag.*l + ~flag.*r;
    rhoWf = op.faceAvg(rhoW);
    rhoGf = op.faceAvg(rhoG);
    
    rhoWf = rhoWf(vtc);
    rhoGf = rhoGf(vtc);
    g = norm(model.gravity);
    t1 = op.connections.faceTopDepth(vtc, 1);
    t2 = op.connections.faceTopDepth(vtc, 2);

    b1 = op.connections.faceBottomDepth(vtc, 1);
    b2 = op.connections.faceBottomDepth(vtc, 2);
    isFine = G.cells.discretization == 1;
    if treatAsCoarse
        % Treat fine cells as VE cells
        if any(isFine)
            % Fine cells: Calculate the water pressure at the TOP [bottom] of the
            % cell by assuming hydrostatic equilibrium within each fine cell.
            % This is then backwards corrected in the "VE" way by further calls
            % to the compute routines.

            T = G.cells.topDepth(isFine);
            B = G.cells.bottomDepth(isFine);
            H = G.cells.height(isFine);
            
            cz = (T + B)/2;
            dz = B - cz;
            %dz = cz - T;
            dpFine = g.*rhoW(isFine).*dz;
            pW(isFine) = pW(isFine) - dpFine;
        end    

        [pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesVE_coarse(model, pW, h, h_max, 1, vtc, rhoW, rhoG, muW, muG, cellsVENot, 'p_entry', opt.p_entry);
        [pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesVE_coarse(model, pW, h, h_max, 2, vtc, rhoW, rhoG, muW, muG, cellsVENot, 'p_entry', opt.p_entry);
      
        isFine_l = false;
        isFine_r = false;
        gdz_g = g*grad(t1, t2);
        gdz_w = g*grad(b1, b2);
    else
        % Treat coarse cells as fine cells
        % Original VE converted to fine cells => change h_max:
        
        [pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesVE_fine(model, pW, h, h_max, 1, vtc, rhoW, rhoG, muW, muG, cellsVENot, 'p_entry', opt.p_entry);       
        [pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesVE_fine(model, pW, h, h_max, 2, vtc, rhoW, rhoG, muW, muG, cellsVENot, 'p_entry', opt.p_entry);

        isFine_l = true;
        isFine_r = true;
        c1 = (t1 + b1)/2;
        c2 = (t2 + b2)/2;
        gdz_g = g*grad(c1, c2);
        gdz_w = gdz_g;
    end
    
    nn = model.operators.N(vtc, :);
    n1 = nn(:,1); n2 = nn(:,2);     

    [pW_l, pG_l, mobW_l, mobG_l] = evaluatePropertiesVE(model, pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l, ...
                                                         isFine_l, n1, h, h_max, vtc, 1, cellsVENot, ...
                                                         've_model', opt.ve_model, 'p_entry', opt.p_entry);
    [pW_r, pG_r, mobW_r, mobG_r] = evaluatePropertiesVE(model, pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r, ...
                                                         isFine_r, n2, h, h_max, vtc, 2, cellsVENot, ...
                                                         've_model', opt.ve_model, 'p_entry', opt.p_entry);    


    dpG   = grad(pG_l, pG_r) - rhoGf .* gdz_g;
    upcg  = (value(dpG)<=0);
    
    vG = -upstr(upcg, mobG_l, mobG_r) .* trans .* dpG;

    dpW   = grad(pW_l, pW_r) - rhoWf .* gdz_w;
    upcw  = (value(dpW)<=0);
    vW = -upstr(upcw, mobW_l, mobW_r) .* trans .* dpW;
end


function [pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, pW_b, sG, h, h_max, H, rhoW, rhoG, muW, muG, ...
                                                      isFine, n_cells, h_global, h_max_global, vtc, index, cellsVENot, varargin)
    % INPUTS:
    %   model: model of type WaterGasMultiVEModel
    %   pW: water pressure
    %   sG: gas saturation
    %   h: thickness of gas plume in coarse/fine/virtual cells
    %   h_max: max reached thickness of gas plume in c/f/v cells
    %   H: thickness of coarse cell
    %   rho[WG]: density for water/gas phase
    %   mu[WG]: mobility for water/gas phase
    %   isFine: boolean indicating if block is fine cell or VE cell
    %   n_cells: coarse cells to evaluate properties for
    %   h_global: global depth of mobile CO2 plume in original coarse cells (i.e. NOT
    %   accounting for virtual cells)
    %   h_max_global: global depth of max CO2 depth reached through history
    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0, 'pc_vetofine', 'fluid');
    opt = merge_options(opt, varargin{:});
    
    g = norm(model.gravity);
    f = model.fluid;
    
    sW = 1 - sG;
    swr = f.krPts.w(1);
    snr = f.krPts.g(1);
    isVE = ~isFine;   
    if numel(isVE) == 1
        isVE = repmat(isVE, numel(n_cells), 1);
    end
      
    n_sealing = ismember(n_cells, model.G.sealingCells);   
    isOriginalVE = model.G.cells.discretization(n_cells) ~= 1; % if coarse cell is originally a VE column (before treating transition regions)    
    
    if isfield(f, 'pcWG') % Fine: brooks-corey. VE: entry pressure (rest). Sealing: fine brooks-corey by default
        pcWG = f.pcWG(sG, n_cells); % should work for both face and cell constraints        
    else
        pcWG = 0;    
    end   
    
    % --- TRY TO CHANGE THIS FOR NVE CELLS TO BE h_max INSTEAD OF h ---
    pcWG_U = (h.*(rhoW - rhoG) - H.*rhoW).*g + opt.p_entry;
    
    if strcmp(opt.ve_model, 'sharp_interface')            
        pcWG = isFine.*pcWG + isVE.*pcWG_U; % old pcWG is p_entry by default for VE cells
    end
    
    if ~isempty(vtc) && ~isempty(index) % transition regions provided - handle them          
       
       t = model.G.cells.topDepth(n_cells);    
       T = model.operators.connections.faceTopDepth(vtc, index);
       b = model.G.cells.bottomDepth(n_cells);
       B = model.operators.connections.faceBottomDepth(vtc, index);       
    
       if all(vtc == (model.operators.connections.veToFineConn | ...
                    model.operators.connections.veTransitionVerticalConn)) && ...
                    all(isVE ~= isOriginalVE)
            % neighbors that are VE cells need to have reconstructed pc, fine
            % cells remain as they are (not included here since isVE == isOriginalVE).            
            z = (T + B)/2; % midpoint of each virtual fine cell
            hp = t + h_global(n_cells); % global depth of plume            
            hmax = t + h_max_global(n_cells);
            %anyCO2 = h_max_global(n_cells) > 0; % no CO2 => entry pressure not reached => pwWG = 0
            
            if strcmp(opt.ve_model, 'sharp_interface')
                % NB: This analytical formula only valid for sharp
                % interface.                
                if strcmp(opt.pc_vetofine, 'satweight')                
%                     pcWG = pG_t + (n_mobile + w_mobile) ...
%                                 + (n_res + w_res) ...
%                                 + (n_brine + w_brine) ...
%                                 - pW;
                    sw_above = (swr.*(z<=hp) + (1-snr).*(z>hp & z<hmax)).*(z-hp); % above largest depth of plume
                    sn_above = ((1-swr).*(z<=hp) + snr.*(z>hp & z<hmax)).*(z-hp);
                    sw_below = ((1-snr).*(hmax-hp) + (z-hmax)).*(z>=hmax); % below largest depth of plume
                    sn_below = snr.*(hmax-hp).*(z>=hmax);
                    
                    pcWG = opt.p_entry - rhoW.*g.*(sw_above + sw_below) ...
                                        - rhoG.*g.*(sn_above + sn_below);
                                    
                elseif strcmp(opt.pc_vetofine, 'hydrostatic')                    
                    pcWG = opt.p_entry - (rhoW-rhoG).*g.*(z-hp);
                    
                elseif strcmp(opt.pc_vetofine, 'new_hydrostatic')
                    pcWG = pcWG_U + rhoG.*g.*z + rhoW.*g.*(H-z);
                end               
            end     
            % OR USE CAPILLARY PRESSURE FUNC ON RECONSTRUCTED SAT ??            
       end
      
    end
    % Gas pressure  
    pG = pW + pcWG;
 
    % Mobility
    swr = f.krPts.w(1);
    snr = f.krPts.g(1);
    %sGmax = h_max./H .* (1-swr); % inverse formula - assumes sharp interface!
    if isprop(model, 'EOSModel')
        krw = f.krO(sW);
    else
        krw = f.krW(sW);
    end
    
    SnMax = (1-swr).*(h_max./H); % h_max local for each cell   
    krg = f.krG(sG, SnMax);         
    
    krwVE = (f.krW(1-snr).*(h_max-h) + (H-h_max))./H; % water mobile in residual zone and pure brine zone
    mobW = (isVE.*krwVE + isFine.*krw)./muW;
    krgVE = f.krG(1-swr).*(h./H); % max CO2 sat in cells of mobile plume is 1-swr (sharp interface assumption)
    mobG = (isVE.*krgVE + isFine.*krg)./muG;
    
    % --- Gas hysteresis ---    
    
    krgVE_PD = @(Sn) f.krG(1-swr, 1-swr) .* Sn./(1-swr); % if Sn=1-swr, then SnMax must also equal 1-swr
    krgVE_PI = @(Sn) f.krG(1-swr, 1-swr) .* max((Sn-snr)./(1-swr-snr), 0); % coarse relperm imbibition
    
    krgVE = @(Sn, SnMax) Hysteresis.Killough(Sn, SnMax, 1-swr, snr, krgVE_PD, krgVE_PI);
                     
    %Sn = (1-swr).*(h./H); % linear function of gas saturation
    Sn = (1-swr-snr).*(h./H) + SnMax.*(snr./(1-swr)); % inverse formula to get coarse saturation
    %mobG = (isVE.*krgVE(Sn, SnMax) + isFine.*krg)./muG; % hysteresis for coarse saturation (indirectly through h, hmax)
    % ----------------------
end

function [vW, vG, mobW, mobG, upcw, upcg] = computeHybridFluxes_test(model, pW, sG, muW, muG, rhoW, rhoG, trans, ...
                                                                            sgMax, varargin)
% Internal function - computes interior fluxes for the hybrid VE models,
% with special treatment for relaxed VE columns.
% 
% PARAMETERS:
%   model   - hybrid model (generated from WaterGasMultiVEModel)
%   pW      - water pressure for hybrid model
%   sG      - co2 saturation (for hybrid model)
%   muW     - water viscosity (cP)
%   muG     - co2 viscosity (cP)
%   rhoW    - water density (kg/m^3)
%   rhoG    - co2 density (kg/m^3)
%   trans   - coarse transmissibilities
%   sgMax   - max co2 saturation reached
%   vG      - accumulated co2 fluxes for each coarse connection
%   vG_smax - accumulated co2 fluxes up to time step when sgMax was reached
%   cB      - all cells connected to top interfaces of semi-perm layers
%   veB     - interfaces (indices) connected to cB
%   cH      - all cells connected to endpoints of semi-perm layers%            
%   veH     - interfaces (indices) connected to cH
%
% RETURNS:
%   vW          - computed water fluxes
%   vG          - computed co2 fluxes
%   mobW        - water mobilites
%   mobG        - co2 mobilites
%   upcw        - upstream flag/bool for water flux
%   upcg        - upstream flag/bool for co2 flux
%   h_global    - (global) height of mobile plume
%   h_T_global  - max height reached for plume residing on top
%   h_B_global  - net height of residual plumes emanating from semi-perm
%                   layers
%   cellsBH     - cell array of cells connected to semi-perm layers
%                   satisfying relaxed VE columns
%               - cellsBH{1}: cB cells, not including cells also cH
%               - cellsBH{2}: cH cells, not including cells also cB
%               - cellsBH{3}: cBH cells, i.e. cells that are both cB and cH
%               - cellsBH{4}: connections associated with cB
%               - cellsBH{5}: connections associated with cH
%               - cellsBH{6}: connections associated with cBH
%   hHi_state   - cell array of heights for cH cells
%               - hHi_state{1}: height of late migrating plume
%               - hHi_state{2}: height of associated virtual VE cell
%   hBHi_state  - cell array of heights for cBH cells
%                   (same elements as hHi_state)
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
    opt = struct('ve_model', 'sharp_interface', 'r', 0); % defaults: sharp interface, negligable horizontal flux for NVEHorz cells (r=0) 
    opt = merge_options(opt, varargin{:});
    
    cellsNVE = [];
    cellsNVEMob = [];
    
    op = model.operators;
    G = model.G;    
    f = model.fluid;
    g = -norm(model.gravity); % !!!
    pv = poreVolume(G, model.rock);
    
    swr = model.fluid.krPts.w(3,2);
    snr = model.fluid.krPts.g(3,2);
    isFine = G.cells.discretization == 1;  
    all_coarse_cells = (1:G.cells.num)';
    
    sgMax_dummy = sgMax.*(snr > 0) + sG.*(snr == 0); % used to assign correct h_max to VE cells filled from bottom    
%     if ~isempty(cellsVE) % once residual plume reaches top, set sgMax to maximum and use standard formulas       
%         sgMax_dummy(cellsVE) = 1-swr; % would otherwise be discontinuous transition from H to snr*H
%     end 
            
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
       
    p_entry = f.pe_rest;    
    
    %[pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h, H, rhoW, rhoG, muW, muG, isFine, isFine);
    [pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h_global, h_max_global, H, rhoW, rhoG, muW, muG, ...
                                                isFine, all_coarse_cells, h_global, h_max_global, {}, {}, ...
                                                've_model', opt.ve_model, 'p_entry', p_entry);
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
    %gdz_g(vIc) = g*(G.cells.topDepth(n1) - G.cells.topDepth(n2));
    % if water pressure computed at bottom of VE cells
    gdz_w(vIc) = g*(G.cells.bottomDepth(n2) - G.cells.bottomDepth(n1));        
    %gdz_w(vIc) = g*(G.cells.bottomDepth(n1) - G.cells.bottomDepth(n2));
    % if water pressure computed at top of VE cells, uncomment:
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
                                                                              transition_ve, true, 'p_entry', p_entry);
    end
    % Treat transition between ve zones in vertical direction (for diffuse
    % leakage)
    transition_vertical = model.operators.connections.veToFineConn | ...
                          model.operators.connections.veTransitionVerticalConn & op.T > 0;
  
    if any(transition_vertical)        
        [vW(transition_vertical), vG(transition_vertical),...
         upcw(transition_vertical), upcg(transition_vertical)] = computeTransitionFluxVE(model, pW, h_global, h_max_global, rhoW, rhoG, muW, muG, ...
                                                                                            transition_vertical, false, 'p_entry', p_entry);
    end

end

function [pW, pW_f, sG, h, h_max, H, rhow, rhog, muw, mug, isFine] = getTransitionValuesVE_coarse(model, pW, h_global, h_max_global, ...
                                                                                                    index, subs, rhoW, rhoG, muW, muG, varargin)    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0, ...
                 'cHorz', [], 'cBottomHorz', [], 'hHi', [], 'hBHi', [], 'Hi', [], 'BHi', []);
    opt = merge_options(opt, varargin{:});
    
    c = model.operators.N(subs, index);
    
    isFine = model.G.cells.discretization(c) == 1;

    t = model.G.cells.topDepth(c);
    T = model.operators.connections.faceTopDepth(subs, index);
    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);

    %sG = getGasSaturationFromHeight(T, t, B, b, h_global(c));   
    swr = model.fluid.krPts.w(3,2);
    snr = model.fluid.krPts.g(3,2);
    % 
    [a_M, a_R, sG] = getGasSatFromHeight_test(T, t, B, b, h_global(c), h_max_global(c), swr, snr);
%     [a_M, a_R, sG] = getGasSatFromHeightFine_test(c, T, t, B, b, h_global(c), h_T_global(c), h_B_global(c), swr, snr, ...
%                                                    opt.cHorz, opt.cBottomHorz, opt.hHi, opt.hBHi, opt.Hi, opt.BHi); % specific for relaxed VE cols
%     
    H = model.operators.connections.faceHeight(subs, index);
    h = a_M.*H;  
    h_max = (a_R + a_M).*H;
    
    g = -norm(model.gravity);
    rhow = rhoW(c);
    rhog = rhoG(c);
    
    pW_f = pW(c);
    muw = muW(c);
    mug = muG(c);      
    
    %pW = getWaterPressureFromHeight(B, t, B, b, h_global(c), h_max_global(c), pW(c), swr, snr, g, rhow, rhog, 'cNVE', cNVE, 'p_entry', opt.p_entry);   
    pW = getPwFromHeight_FF(B, t, b, h_global(c), h_max_global(c), pW(c), g, rhow, rhog, swr, snr, ...
                                'cNVE', [], 'p_entry', opt.p_entry, 'res_type', []);
end

function [pW, pW_f, sG, h, h_max, H, rhow, rhog, muw, mug, isFine] = getTransitionValuesVE_fine(model, pW, h_global, h_max_global, index, subs, rhoW, rhoG, muW, muG, varargin)    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0, ...
                 'cHorz', [], 'cBottomHorz', [], 'hHi', [], 'hBHi', [], 'Hi', [], 'BHi', []);
    opt = merge_options(opt, varargin{:});
    
    c = model.operators.N(subs, index);
    t = model.G.cells.topDepth(c);
    
    isFine = model.G.cells.discretization(c) == 1;

    T = model.operators.connections.faceTopDepth(subs, index);

    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);
    H = model.operators.connections.faceHeight(subs, index);
    
    %sG = getGasSaturationFromHeight(T, t, B, b, h_global(c));
    swr = model.fluid.krPts.w(3,2);
    snr = model.fluid.krPts.g(3,2);    
       
    
    [a_M, a_R, sG] = getGasSatFromHeight_test(T, t, B, b, h_global(c), h_max_global(c), swr, snr);   

    h = a_M.*H;  
    h_max = (a_R + a_M).*H;
    
    g = -norm(model.gravity);
    rhow = rhoW(c);
    rhog = rhoG(c);
    C = (T + B)/2;
    pW_f = pW(c);
    muw = muW(c);
    mug = muG(c);    

    pW_c = getPwFromHeight_FF(C, t, b, h_global(c), h_max_global(c), pW_f, g, rhow, rhog, swr, snr, ...
                            'cNVE', [], 'p_entry', opt.p_entry, 'res_type', []);
    pW = isFine.*pW_f + ~isFine.*pW_c;
end


function [vW, vG, upcw, upcg] = computeTransitionFluxVE(model, pW, h, h_max, rhoW, rhoG, muW, muG, vtc, treatAsCoarse, varargin)
    % INPUTS:
    %   h: global depth of mobile plume
    %   h_max: global height of residual part
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0, ...
                 'cHorz', [], 'cBottomHorz', [], 'hHi', [], 'hBHi', [], 'Hi', [], 'BHi', []);
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
    g = -norm(model.gravity);
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
            dz = cz - B; % B - cz            
            dpFine = g.*rhoW(isFine).*dz;
            pW(isFine) = pW(isFine) - dpFine;
        end    

        [pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesVE_coarse(model, pW, h, h_max, 1, vtc, rhoW, rhoG, muW, muG, 'p_entry', opt.p_entry);
        [pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesVE_coarse(model, pW, h, h_max, 2, vtc, rhoW, rhoG, muW, muG, 'p_entry', opt.p_entry);
      
        isFine_l = false;
        isFine_r = false;
        gdz_g = g*grad(t1, t2);
        gdz_w = g*grad(b1, b2);
    else
        % Treat coarse cells as fine cells
        % Original VE converted to fine cells => change h_max:
        
        [pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesVE_fine(model, pW, h, h_max, 1, vtc, rhoW, rhoG, muW, muG, 'p_entry', opt.p_entry);
        [pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesVE_fine(model, pW, h, h_max, 2, vtc, rhoW, rhoG, muW, muG, 'p_entry', opt.p_entry);

        isFine_l = true;
        isFine_r = true;
        c1 = (t1 + b1)/2;
        c2 = (t2 + b2)/2;
        gdz_g = g*grad(c1, c2);
        gdz_w = gdz_g;
    end
    
    nn = model.operators.N(vtc, :);
    n1 = nn(:,1); n2 = nn(:,2);     

    [pW_l, pG_l, mobW_l, mobG_l] = evaluatePropertiesVE(model, pW_l, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l, ...
                                                         isFine_l, n1, h, h_max, vtc, 1, ...
                                                         've_model', opt.ve_model, 'p_entry', opt.p_entry);
    [pW_r, pG_r, mobW_r, mobG_r] = evaluatePropertiesVE(model, pW_r, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r, ...
                                                         isFine_r, n2, h, h_max, vtc, 2, ...
                                                         've_model', opt.ve_model, 'p_entry', opt.p_entry);    


    dpG   = grad(pG_l, pG_r) - rhoGf .* gdz_g;
    upcg  = (value(dpG)<=0);
    
    vG = -upstr(upcg, mobG_l, mobG_r) .* trans .* dpG;

    dpW   = grad(pW_l, pW_r) - rhoWf .* gdz_w;
    upcw  = (value(dpW)<=0);
    vW = -upstr(upcw, mobW_l, mobW_r) .* trans .* dpW;
end


function [pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h, h_max, H, rhoW, rhoG, muW, muG, ...
                                                      isFine, n_cells, h_global, h_max_global, vtc, index, varargin)
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
    
    %g = -norm(model.gravity);
    g = norm(model.gravity);
    f = model.fluid;
    
    sW = 1 - sG;
    swr = f.krPts.w(3,2);
    snr = f.krPts.g(3,2);
    isVE = ~isFine;   
    if numel(isVE) == 1
        isVE = repmat(isVE, numel(n_cells), 1);
    end
      
    n_sealing = ismember(n_cells, model.G.sealingCells);   
    isOriginalVE = model.G.cells.discretization(n_cells) ~= 1; % if coarse cell is originally a VE column (before treating transition regions)    
    
    if isfield(f, 'pcWG') % Fine: brooks-corey. VE: entry pressure (rest). Sealing: fine brooks-corey by default
        % NB: Change cell-indexation to be correct region (not just 1)!
        pcWG = f.pcWG{1}(sG); % should work for both face and cell constraints        
    else
        pcWG = 0;    
    end   
    
    % --- TRY TO CHANGE THIS FOR NVE CELLS TO BE h_max INSTEAD OF h ---
    pcWG_U = (h.*(rhoW - rhoG) - H.*rhoW).*g + opt.p_entry;
    
    if strcmp(opt.ve_model, 'sharp_interface')            
        pcWG = isFine.*pcWG + isVE.*pcWG_U; % old pcWG is p_entry by default for VE cells
    end
        
    % Gas pressure  
    pG = pW + pcWG;
 
    % Mobility
    if isprop(model, 'EOSModel')
        krw = f.krO(sW);
    else
        krw = f.krW{3}(sW);
        %for fac=unique(model.G.facies)
        %krw(model.G.facies == fac) = f.krw{fac}(sW(model.G.facies == fac));
    end
    
    SnMax = (1-swr).*(h_max./H); % h_max local for each cell   
    %krg = f.krG(sG, SnMax);      
    krg = f.krG{3}(sG);
    
    krwVE = (f.krW{3}(1-snr).*(h_max-h) + (H-h_max))./H; % water mobile in residual zone and pure brine zone
    mobW = (isVE.*krwVE + isFine.*krw)./muW;
    krgVE = f.krG{3}(1-swr).*(h./H); % max CO2 sat in cells of mobile plume is 1-swr (sharp interface assumption) 
    mobG = (isVE.*krgVE + isFine.*krg)./muG;
    
end


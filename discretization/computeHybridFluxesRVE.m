function [vW, vG, mobW, mobG, upcw, upcg, ...
            h_global, h_max_global, h_T_global, h_B_global, ...
            cellsBH, cHorz_state, cBottomHorz_state] = computeHybridFluxesRVE(model, pW, sG, muW, muG, rhoW, rhoG, trans, ...
                                                                            sgMax, vG, vG_smax, cB, veB, cH, veH, varargin)
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
%   cHorz_state   - cell array of heights for cH cells
%                 - cHorz_state{1}: height of late migrating plume
%                 - cHorz_state{2}: height of associated virtual VE cell
%   cBottomHorz_state  - cell array of heights for cBH cells
%                       (same elements as cHorz_state)
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
    
    op = model.operators;
    G = model.G;    
    f = model.fluid;
    g = norm(model.gravity);
    pv = poreVolume(G, model.rock);
    
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
    isFine = G.cells.discretization == 1;  
    all_coarse_cells = (1:G.cells.num)';
    
    sgMax_dummy = sgMax.*(snr > 0) + sG.*(snr == 0); % used to assign correct h_max to VE cells filled from bottom    
            
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

    % --------------------------   
    h_T_global = h_max_global.*(snr > 0) + h_global.*(snr == 0);
    h_B_global = H;    
               
    % --- Special treatment of RVE cells with bottom/horizontal fluxes from sealing layers ---   
    [h_ve, hmax_ve, hT_ve, hB_ve, ...
        cellsBH, cHorz_state, cBottomHorz_state] = getRVEHeights(model, cB, veB, cH, veH, sG, sgMax, vG, vG_smax, H, swr, snr, pv);
    
    cBottom = cellsBH{1};
    cHorz = cellsBH{2};
    cBottomHorz = cellsBH{3};
    cB = unique(cBottom);
    cH = unique(cHorz);
    cBH = unique(cBottomHorz);
    
    % ----------------------------
    h_global(cB) = h_ve.B;
    h_max_global(cB) = hmax_ve.B;
    h_T_global(cB) = hT_ve.B;
    h_B_global(cB) = hB_ve.B;
    
    h_global(cH) = h_ve.H;
    h_max_global(cH) = hmax_ve.H;
    h_T_global(cH) = hT_ve.H;
    h_B_global(cH) = hB_ve.H;
    
    h_global(cBH) = h_ve.BH; 
    h_max_global(cBH) = hmax_ve.BH;
    h_T_global(cBH) = hT_ve.BH;
    h_B_global(cBH) = hB_ve.BH;
    % ----------------------------   
       
    p_entry = f.pe_rest;    

    [pW, pG, mobW, mobG] = evaluatePropertiesRVE(model, pW, sG, h_global, h_max_global, H, rhoW, rhoG, muW, muG, ...
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
    gdz_w(vIc) = g*(G.cells.bottomDepth(n2) - G.cells.bottomDepth(n1)); 
    
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
         upcw(transition_ve), upcg(transition_ve)] = computeTransitionFluxRVE(model, pW, h_global, h_max_global, h_T_global, h_B_global, rhoW, rhoG, muW, muG, ...
                                                                              transition_ve, true, 'p_entry', p_entry, ...
                                                                              'cHorz', cHorz, 'cBottomHorz', cBottomHorz, ...
                                                                              'hHi', cHorz_state{1}, 'hBHi', cBottomHorz_state{1}, ...
                                                                              'Hi', cHorz_state{2}, 'BHi', cBottomHorz_state{2});
    end
    % Treat transition between ve zones in vertical direction (for diffuse
    % leakage)
    transition_vertical = model.operators.connections.veToFineConn | ...
                          model.operators.connections.veTransitionVerticalConn & op.T > 0;
  
    if any(transition_vertical)        
        [vW(transition_vertical), vG(transition_vertical),...
         upcw(transition_vertical), upcg(transition_vertical)] = computeTransitionFluxRVE(model, pW, h_global, h_max_global, h_T_global, h_B_global, rhoW, rhoG, muW, muG, ...
                                                                                            transition_vertical, false, 'p_entry', p_entry, ...
                                                                                            'cHorz', cHorz, 'cBottomHorz', cBottomHorz, ...
                                                                                            'hHi', cHorz_state{1}, 'hBHi', cBottomHorz_state{1}, ...
                                                                                            'Hi', cHorz_state{2}, 'BHi', cBottomHorz_state{2});
    end

end

function [pW, pW_f, sG, h, h_max, H, rhow, rhog, muw, mug, isFine] = getTransitionValuesRVE_coarse(model, pW, h_global, h_max_global, h_T_global, h_B_global, ...
                                                                                                    index, subs, rhoW, rhoG, muW, muG, varargin)    
    % Compute values at specified virtual cells for a (coarse) VE-to-VE
    % horizontal transition zone.
    % INPUTS:
    %   model: hybrid model (e.g., instance of WaterGasHybridModel)
    %   pW: water pressure
    %   h_global: height of mobile plumes in each cell of hybrid grid
    %   h_max_global: maximum reached height of mobile plumes
    %   h_T_global: maximum reached height of buoyantly segragated plume
    %   h_B_global: net height of plumes originating from sealing layers
    %   index: the neighboring cell to get values from
    %   subs: subset of global cells part of transition zone
    %   rhoX: density of phase X
    %   muX: viscosity of phase X
    % RETURNS:
    %   pW: water pressure evaluated for virtual cells of subset
    %   pW_f: water pressure evaluated for fine cells of subset
    %   sG: gas saturation
    %   h: height of mobile plume
    %   h_max: net height of residual parts of plumes
    %   H: height of cells in subset
    %   rhoX: density of phase X
    %   muX: viscosity of phase X
    %   isFine: true if cells is fine, false if VE

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
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
    % 
    [a_M, a_R, sG] = getGasSatFromHeightVirtual(T, t, B, b, h_global(c), h_T_global(c), h_B_global(c), swr, snr);
  
    H = model.operators.connections.faceHeight(subs, index);
    h = a_M.*H;  
    h_max = (a_R + a_M).*H;
    
    g = norm(model.gravity);
    rhow = rhoW(c);
    rhog = rhoG(c);
    
    pW_f = pW(c);
    muw = muW(c);
    mug = muG(c);      
    
    pW = getPwFromHeight(B, t, b, h_global(c), h_T_global(c), pW(c), g, rhow, rhog, swr, snr, ...
                                'cNVE', [], 'p_entry', opt.p_entry(c), 'res_type', []); % opt.p_entry
end

function [pW, pW_f, sG, h, h_max, H, rhow, rhog, muw, mug, isFine] = getTransitionValuesRVE_fine(model, pW, h_global, h_max_global, h_T_global, h_B_global, ...
                                                                                                    index, subs, rhoW, rhoG, muW, muG, varargin)    
    % Compute values at specified virtual cells for a (fine) VE-to-VE
    % vertical transition zone or VE-to-fine transition zone.
    % INPUTS and OUTPUTS:
    %   see getTransitionValuesRVE_coarse
    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0, ...
                 'cHorz', [], 'cBottomHorz', [], 'hHi', [], 'hBHi', [], 'Hi', [], 'BHi', []);
    opt = merge_options(opt, varargin{:});
    
    c = model.operators.N(subs, index);   
    
    isFine = model.G.cells.discretization(c) == 1;

    t = model.G.cells.topDepth(c);
    T = model.operators.connections.faceTopDepth(subs, index);
    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);
    H = model.operators.connections.faceHeight(subs, index);
    
    %sG = getGasSaturationFromHeight(T, t, B, b, h_global(c));
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);    
  
    [a_M, a_R, sG] = getGasSatFromHeightVirtual(T, t, B, b, h_global(c), h_T_global(c), h_B_global(c), swr, snr);   

    h = a_M.*H;  
    h_max = (a_R + a_M).*H;
    
    g = norm(model.gravity);
    rhow = rhoW(c);
    rhog = rhoG(c);
    C = (T + B)/2;
    pW_f = pW(c);
    muw = muW(c);
    mug = muG(c);    
   
    pW_c = getPwFromHeight(C, t, b, h_global(c), h_T_global(c), pW_f, g, rhow, rhog, swr, snr, ...
                            'cNVE', [], 'p_entry', opt.p_entry(c), 'res_type', []); % opt.p_entry -> opt.p_entry(c)
    pW = isFine.*pW_f + ~isFine.*pW_c;
end


function [vW, vG, upcw, upcg] = computeTransitionFluxRVE(model, pW, h, h_max, h_T, h_B, rhoW, rhoG, muW, muG, vtc, treatAsCoarse, varargin)
    % Compute discrete flux at interface shared by two cells where either
    % one or both are RVE.
    % INPUTS:
    %   model: hybrid model
    %   pW: water pressure
    %   h: height of mobile plume
    %   h_max: net height of residual plumes
    %   h_T: maximum height reached of top segregated plume
    %   h_B: net height reached by isolated plumes from sealing layers
    %   rhoX: density of phase X
    %   muX: viscosity of phase X
    %   vtc: connections/faces to compute flux for
    %   treatAsCoarse: true if virtual cells are treated as VE cells
    % RETURNS:
    %   vW: flux of water phase
    %   vG: flux of gas phase
    %   upcw: upwind operators for water phase
    %   upcg: upwind operator for gas phase

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

        [pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesRVE_coarse(model, pW, h, h_max, h_T, h_B, 1, vtc, rhoW, rhoG, muW, muG, 'p_entry', opt.p_entry, ...
                                                                                                            'cHorz', opt.cHorz, 'cBottomHorz', opt.cBottomHorz, ...
                                                                                                            'hHi', opt.hHi, 'hBHi', opt.hBHi, ...
                                                                                                            'Hi', opt.Hi, 'BHi', opt.BHi);
        [pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesRVE_coarse(model, pW, h, h_max, h_T, h_B, 2, vtc, rhoW, rhoG, muW, muG, 'p_entry', opt.p_entry, ...
                                                                                                            'cHorz', opt.cHorz, 'cBottomHorz', opt.cBottomHorz, ...
                                                                                                            'hHi', opt.hHi, 'hBHi', opt.hBHi, ...
                                                                                                            'Hi', opt.Hi, 'BHi', opt.BHi);
      
        isFine_l = false;
        isFine_r = false;
        gdz_g = g*grad(t1, t2);
        gdz_w = g*grad(b1, b2);
    else
        % Treat coarse cells as fine cells
        % Original VE converted to fine cells => change h_max:
        
        [pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesRVE_fine(model, pW, h, h_max, h_T, h_B, 1, vtc, rhoW, rhoG, muW, muG, 'p_entry', opt.p_entry, ...
                                                                                                           'cHorz', opt.cHorz, 'cBottomHorz', opt.cBottomHorz, ...
                                                                                                            'hHi', opt.hHi, 'hBHi', opt.hBHi, ...
                                                                                                            'Hi', opt.Hi, 'BHi', opt.BHi);
        [pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesRVE_fine(model, pW, h, h_max, h_T, h_B, 2, vtc, rhoW, rhoG, muW, muG, 'p_entry', opt.p_entry, ...
                                                                                                            'cHorz', opt.cHorz, 'cBottomHorz', opt.cBottomHorz, ...
                                                                                                            'hHi', opt.hHi, 'hBHi', opt.hBHi, ...
                                                                                                            'Hi', opt.Hi, 'BHi', opt.BHi);

        isFine_l = true;
        isFine_r = true;
        c1 = (t1 + b1)/2;
        c2 = (t2 + b2)/2;
        gdz_g = g*grad(c1, c2);
        gdz_w = gdz_g;
    end
    
    nn = model.operators.N(vtc, :);
    n1 = nn(:,1); n2 = nn(:,2);     

    [pW_l, pG_l, mobW_l, mobG_l] = evaluatePropertiesRVE(model, pW_l, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l, ...
                                                         isFine_l, n1, h, h_max, vtc, 1, ...
                                                         've_model', opt.ve_model, 'p_entry', opt.p_entry);
    [pW_r, pG_r, mobW_r, mobG_r] = evaluatePropertiesRVE(model, pW_r, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r, ...
                                                         isFine_r, n2, h, h_max, vtc, 2, ...
                                                         've_model', opt.ve_model, 'p_entry', opt.p_entry);    


    dpG   = grad(pG_l, pG_r) - rhoGf .* gdz_g;
    upcg  = (value(dpG)<=0);
    
    vG = -upstr(upcg, mobG_l, mobG_r) .* trans .* dpG;

    dpW   = grad(pW_l, pW_r) - rhoWf .* gdz_w;
    upcw  = (value(dpW)<=0);
    vW = -upstr(upcw, mobW_l, mobW_r) .* trans .* dpW;
end


function [pW, pG, mobW, mobG] = evaluatePropertiesRVE(model, pW, sG, h, h_max, H, rhoW, rhoG, muW, muG, ...
                                                      isFine, n_cells, h_global, h_max_global, vtc, index, varargin)
    % Evaluate properties for virtual cells in RVE region.
    % Used to compute fluxes at associated interfaces.
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
    %   n_cells: cells (from hybrid model) to evaluate properties for
    %   h_global: global depth of mobile CO2 plume in original coarse cells (i.e. NOT
    %   accounting for virtual cells)
    %   h_max_global: global depth of max CO2 depth reached through history
    %   vtc: connections/faces to evaluate properties for
    %   index: index of neighboring cells in subset
    % RETURNS:
    %   pW: water pressure
    %   pG: gas pressure
    %   mobW: water mobility
    %   mobG: gas mobility
    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0, 'pc_vetofine', 'fluid');
    opt = merge_options(opt, varargin{:});
    
    g = norm(model.gravity);
    f = model.fluid;
    rock = model.rock;
    
    sW = 1 - sG;
    isVE = ~isFine;   
    if numel(isVE) == 1
        isVE = repmat(isVE, numel(n_cells), 1);
    end
          
    if isfield(f, 'pcWG') % Fine: brooks-corey. VE: entry pressure (rest). Sealing: fine brooks-corey by default
        pcWG = f.pcWG(sG, n_cells); % should work for both face and cell constraints        
    else
        pcWG = 0;    
    end   
    
    pcWG_U = (h.*(rhoW - rhoG) - H.*rhoW).*g + opt.p_entry(n_cells); % + opt.p_entry 
    
    if strcmp(opt.ve_model, 'sharp_interface')            
        pcWG = isFine.*pcWG + isVE.*pcWG_U; % old pcWG is p_entry by default for VE cells
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
            
    krg = f.krG(sG);
    
    krwVE = (f.krW(1-snr).*(h_max-h) + (H-h_max))./H; % water mobile in residual zone and pure brine zone
    mobW = (isVE.*krwVE + isFine.*krw)./muW;
    krgVE = f.krG(1-swr).*(h./H); % max CO2 sat in cells of mobile plume is 1-swr (sharp interface assumption) 
    mobG = (isVE.*krgVE + isFine.*krg)./muG;
    
end

function [h, h_max, h_T, h_B, cellsBH, ...
            h_Hi_state, h_BHi_state] = getRVEHeights(model, cB_all, veB_all, cH_all, veH_all, sG, sgMax, vG, vG_smax, H, swr, snr, pv)
    % Compute plume heights for relaxed VE columns.
    %
    % PARAMETERS:
    %   model: instance of WaterGasHybridModel
    %   cB_all: RVE columns above bottom sealing layer
    %   veB_all: connections between cB_all and associated sealing layer
    %   cH_all: RVE columns laterally adjacent to a sealing layer
    %   veH_all: connections between cH_all and associated sealing layer
    %   sG: gas saturation
    %   sgMax: max gas saturation reached
    %   vG: gas phase flux
    %   vG_smax: gas phase flux at time step where sgMax was reached
    %   H: height of cells in hybrid model
    %   swr: water residual saturation
    %   snr: gas residual saturation
    %   pv: pore volume
    %
    % RETURNS:
    %   h: struct with heights of mobile CO2 in RVE cells
    %   h_max: struct with *net* maximum height reached for mobile CO2 in
    %           RVE cells (including top buoyant plume and late-stage plumes)
    %   h_T: height of maximum height reached for buoyantly segregated
    %           plume on top of RVE cells
    %   h_B: *net* height of late-stage migrating plumes in RVE cells
    %   cellsBH: cell array of RVE cells and associated connections,
    %           distinguished into 1) only transition to bottom sealing layers, 2)
    %           only transitions to horizontal sealing layers, and 3) transitions
    %           to bottom AND horizontal sealing layers.
    %   h_Hi_state: cell array of plume heights (hi) and cell heights (Hi) for RVE
    %               columns represented by horizontal transition to 
    %               sealing layers. 
    %               Needed for fine-scale reconstruction.
    %   h_BHi_state: cell array of plume heights (hi) and cell heights (BHi) for RVE
    %               columns represented by horizontal AND bottom transitions to 
    %               sealing layers.
    %               Needed for fine-scale reconstruction.
    %
    %   Height variables h, h_max, h_T and h_B have the following fields:
    %       h*.B: heights for RVE columns transitioning to BOTTOM sealing
    %             layer.
    %       h*.H: heights for RVE columns transitioning to HORIZONTAL
    %             sealing layer(s).
    %       h*.BH: heights for RVE columns transitioning to HORIZONTAL
    %              *and* BOTTOM sealing layers.

    op = model.operators;
    p = model.G.partition;
    veBottom = op.connections.veToFineVertical | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;
    
    veHorz = model.operators.connections.veTransitionHorizontalConn;    
    
    cH_all = [cH_all{1}; cH_all{2}]; % merge neighboring cells of veHorz transition interface
    veH_all = [veH_all{1}; veH_all{2}]; % merge connection indices corresponding to neighboring cells of veHorz
    
    h = struct;
    h_max = struct;
    h_T = struct;
    h_B = struct;
    
    cellsBH = {};
    
    % --- BOTTOM ---
    cB = ~ismember(cB_all, cH_all); % VE bottom fluxes but no veHorizontal fluxes
    veB = veB_all(cB);
    cB = cB_all(cB);
    cellsBH = cat(1, cellsBH, {cB, veB});
    
    vG_bottom = vG(veBottom); % take abs val since fluxes are negative from bottom and up
    vG_bottom = abs(vG_bottom(veB));
    vG_bottom_smax = vG_smax(veBottom); % global indexing
    vG_bottom_smax = abs(vG_bottom_smax(veB)); % local indexing
    
    S = sG(cB); % VE bottom cells cB used directly as index because sG is ordered by cell numbers by default (1:G.cells.num)
    sMax = sgMax(cB).*(snr > 0) + S.*(snr == 0);
    
    S_B = 1./pv(cB) .* vG_bottom; % CO2 sat originating from bottom flux
    S_B_smax = 1./pv(cB) .* vG_bottom_smax; % earlier CO2 sat from bottom flux at time where smax was reached --> gives correct max saturation for top co2 part
    S_max_T = sMax - S_B_smax.*(snr > 0); % max saturation reached for top part (excluding part originating from bottom flux)
    
    h_T.B = max(H(cB).*S_max_T./(1-swr), 0); % max to avoid negative discrepancies
    % For h_B we use current saturation at bottom S_B, since
    % this is independent of smax
    h_B.B = max(h_T.B, H(cB).*(1-(S_B./(snr+eps)).*(snr>0))); % h_B will become h_T before h_T becomes negative (S_B > sMax), then h_max_global = H and col is in VE                 
    
    Snr_tot = (snr./H(cB)).*(h_T.B + (H(cB) - h_B.B)); % total residual saturation in column (NB: includes residual part in mobile zone!)
    S_mob = S - Snr_tot; % mobile saturation (at top of column)
       
    h.B = H(cB).*(S_mob./(1-swr-snr)); % subtract snr since this was included in Snr_tot
    h_max.B = H(cB) - h_B.B + h_T.B;    
    
    
    
    % --- HORIZONTAL ---
    cH = ~ismember(cH_all, cB_all);
    veH = veH_all(cH);
    cH = cH_all(cH);
    cellsBH = cat(1, cellsBH, {cH, veH});
    
    vG_horz = vG(veHorz); % same size as veH_all{1}/veH_all{2} 
    vG_horz = abs(vG_horz(veH));   
    vG_horz_smax = vG_smax(veHorz);
    vG_horz_smax = abs(vG_horz_smax(veH));
    
    cH_u = unique(cH);        
       
    t = model.G.cells.topDepth(cH);
    T1 = model.operators.connections.faceTopDepth(:, 1);
    T2 = model.operators.connections.faceTopDepth(:, 2);
    T = [T1(veHorz); T2(veHorz)]; % vertically concatenated    
       
    T = T(veH);
    H_i = T - t; % height from top of virtual ve cell to top of parent ve cell        
         
    % --- Bottom saturations for each virtual cell ---
    S = sG(cH); 
    sMax = sgMax(cH).*(snr > 0) + S.*(snr == 0);   
    S_B = 1./pv(cH) .* vG_horz; % CO2 sat originating from bottom flux  
    S_B_smax = 1./pv(cH) .* vG_horz_smax; % earlier CO2 sat from bottom flux at time where smax was reached --> gives correct max saturation for top co2 part
    S_B_smax = accumarray(cH, S_B_smax);
    S_B_smax_acc = S_B_smax(cH); % choose all cells
    
    S_max_T = sMax - S_B_smax_acc.*(snr > 0); % will be the same for same cells in cH_all
    
    S_Bi = S_B./((H_i./H(cH).*(snr+eps))); % scaled bottom saturation for VE virtual cells
    h_Bi = min(H_i.*(1-S_Bi.*(snr>0)), H_i); % height from top of virtual cell to bottom plume
    h_Ti = max(H(cH).*S_max_T./(1-swr), 0); % height from top of virtual cell to top plume
    
    % Accumulate bottom saturations for virtual cells of same VE col       
    S = sG(cH_u);
    sMax = sgMax(cH_u).*(snr > 0) + S.*(snr == 0); 
    S_B_smax_acc = S_B_smax(cH_u); % only choose unique cells
    S_max_T = sMax - S_B_smax_acc.*(snr > 0); % max saturation reached for top part (excluding part originating from bottom flux)
    
    h_T.H = max(H(cH_u).*S_max_T./(1-swr), 0); % max to avoid negative discrepancies
    h_Bi = max(h_Ti, h_Bi); % h_Ti or h_T.H ???
    
    h_Hi_state = {h_Bi, H_i};
    
    % Accumulate bottom heights  
    dh_Bi = accumarray(cH, value(H_i - h_Bi)); % THICKNESS of bottom plume in virtual cell i  
    h_B.H = max(h_T.H, H(cH_u) - dh_Bi(cH_u)); % h_B will become h_T before h_T becomes negative (S_B > sMax), then h_max_global = H and col is in VE                 
    % -------------------------------------
    
    Snr_tot = (snr./H(cH_u)).*(h_T.H + (H(cH_u) - h_B.H)); % total residual saturation in column (NB: includes residual part in mobile zone!)
    S_mob = S - Snr_tot; % mobile saturation (at top of column)
       
    % mobile and max height function of accumulated bottom sat
    h.H = H(cH_u).*(S_mob./(1-swr-snr)); % subtract snr since this was included in Snr_tot
    h_max.H = H(cH_u) - h_B.H + h_T.H;
         
    
    % --- BOTTOM_HORIZONTAL ---
    cBH_bottom = ismember(cB_all, cH_all);   
    veB = veB_all(cBH_bottom); % ???
    cBH_horz = ismember(cH_all, cB_all);
    veH_all = veH_all(cBH_horz);
    
    cBH = cB_all(cBH_bottom); % desired cells from bottom list. NB: only one bottom connection per cell!
    cBH_h = cH_all(cBH_horz); % desired cells from horizontal list (may contain duplicates if multiple horz connections to a given cell). doesnt matter if we choose to index by cH_all or cB (but if choosing cB, bottom_and_horz need to be calculated for veH_all first and veB last)   
    [~, ~, cBH_h_idx] = unique(cBH_h);
    cellsBH = cat(1, cellsBH, {[cBH_h; cBH], [veH_all; veB]});
    
    % Now do same calculations, but sum bottom heights from bottom fluxes
    % and horizontal fluxes into a coarse h_B and THEN calculate h_T, h and
    % h_max.      
       
    vG_bottom = vG(veBottom); % take abs val since fluxes are negative from bottom and up
    vG_bottom = abs(vG_bottom(veB));
    vG_bottom_smax = vG_smax(veBottom); % global indexing
    vG_bottom_smax = abs(vG_bottom_smax(veB)); % local indexing
       
    vG_horz = vG(veHorz); % same size as veH_all{1}/veH_all{2} 
    vG_horz = abs(vG_horz(veH_all));   
    vG_horz_smax = vG_smax(veHorz);
    vG_horz_smax = abs(vG_horz_smax(veH_all));

    % --- Bottom saturations for horizontal VE flux ---    
    S = sG(cBH_h); 
    sMax = sgMax(cBH_h).*(snr > 0) + S.*(snr == 0);
    
    t = model.G.cells.topDepth(cBH_h);
    T = [model.operators.connections.faceTopDepth(veHorz, 1);
         model.operators.connections.faceTopDepth(veHorz, 2)];
    T = T(veH_all);
    H_i = T - t; % height from top of virtual ve cell to top of parent ve cell                         
    
    S_B = 1./pv(cBH_h) .* vG_horz; % CO2 sat originating from bottom flux    
    S_Bh_smax = 1./pv(cBH_h) .* vG_horz_smax; % earlier CO2 sat from bottom flux at time where smax was reached --> gives correct max saturation for top co2 part
    S_Bh_smax = accumarray(cBH_h, S_Bh_smax);
    S_Bh_smax_acc = S_Bh_smax(cBH_h); % choose all cells (duplicates allowed)
    S_max_Th = sMax - S_Bh_smax_acc.*(snr > 0);
    
    S_Bh = S_B./(H_i./H(cBH_h).*(snr+eps)); % scaled bottom saturation for VE virtual cells
    h_Th = max(H(cBH_h).*S_max_Th./(1-swr), 0); % depth of top plume for each virtual VE cell
    h_Bh = min(H_i.*(1-S_Bh.*(snr>0)), H_i); % depth of bottom plume for each virtual VE cell   
    %h_Bh = max(h_Th, h_Bh);    
    
    % --- Bottom saturation from bottom flux --- 
    S = sG(cBH);
    sMax = sgMax(cBH).*(snr > 0) + S.*(snr == 0); 
    
    S_B = 1./pv(cBH) .* vG_bottom; % CO2 sat originating from bottom flux
    S_Bb_smax = 1./pv(cBH) .* vG_bottom_smax; % earlier CO2 sat from bottom flux at time where smax was reached --> gives correct max saturation for top co2 part
    % No need to accumulate for bottom fluxes - only one connection per cell anyway! 
    S_max_Tb = sMax - S_Bb_smax.*(snr > 0); % max saturation reached for top part (excluding part originating from bottom flux)
    
    h_Tb = max(H(cBH).*S_max_Tb./(1-swr), 0); % max to avoid negative discrepancies
    % For h_B we use current saturation at bottom S_B, since
    % this is independent of smax
    h_Bb = min(H(cBH).*(1-(S_B./(snr+eps)).*(snr>0)), H(cBH));    
    %h_Bb = max(h_Tb, h_Bb); % h_B will become h_T before h_T becomes negative (S_B > sMax), then h_max_global = H and col is in VE                                
    % -----------------------------------------------
        
    % --- Accumulate bottom saturations for virtual cells of same VE col,
    % from BOTH horizontal and bottom transitions ---
    S = sG(cBH); % indexing by cBH chooses unique cells satisfying bottom and horz fluxes
    sMax = sgMax(cBH).*(snr > 0) + S.*(snr == 0); 
    
    S_B_smax_acc = S_Bh_smax(cBH) + S_Bb_smax; % add residual plumes from horizontal and bottom transitions 
    
    S_max_T = sMax - S_B_smax_acc.*(snr > 0); % max saturation reached for top part (excluding part originating from bottom flux)    
    h_T.BH = max(H(cBH).*S_max_T./(1-swr), 0); % max to avoid negative discrepancies
    
    % Cap h_Bh and h_Bb to not exceed global depth of top residual plume
    h_Bh = max(h_Bh, h_T.BH(cBH_h_idx)); % duplicate height for repeated cells
    h_Bb = max(h_Bb, h_T.BH);
    
    % Store in state
    h_BHi_state = {[h_Bh; h_Bb], [H_i; H(cBH)]}; % add all residual plume heights (bottom + horz) into one array
    
    % Accumulate residual plume heights     
    dh_Bh = accumarray(cBH_h, value(H_i - h_Bh));
    dh_Bh = dh_Bh(cBH); % choose unique cells cBH, since we want only one value for accumulated residual plume heights per VE cell
    dh_B = dh_Bh + value(H(cBH)-h_Bb);    
    h_B.BH = max(h_T.BH, H(cBH) - dh_B); % h_B will become h_T before h_T becomes negative (S_B > sMax), then h_max_global = H and col is in VE                 
    % -------------------------------------
     
    Snr_tot = (snr./H(cBH)).*(h_T.BH + (H(cBH) - h_B.BH)); % total residual saturation in column (NB: includes residual part in mobile zone!)
    S_mob = S - Snr_tot; % mobile saturation (at top of column)
       
    % mobile and max height function of accumulated bottom sat
    h.BH = H(cBH).*(S_mob./(1-swr-snr)); % subtract snr since this was included in Snr_tot
    h_max.BH = H(cBH) - h_B.BH + h_T.BH;
end
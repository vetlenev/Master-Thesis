function [vW, vG, mobW, mobG, upcw, upcg] = computeHybridFluxesVEres_test(model, pW, sG, muW, muG, rhoW, rhoG, trans, sgMax, vG_B, vG_B_smax, cB, veB, cellsNVE, cellsNVEMob, cellsVE, cellsNVEHorz, varargin)
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
  
%     if ~isempty(cellsNVE) % only change height for transition zones
%        % NB: Coarse saturation remains unchanged, but the DISTRIBUTION
%        % Mass conservation represented by sG, not sgMax!
%        h_global(cellsNVE) = 0; % no mobile part of plume => no CO2 mobility!      
%        h_max_global(cellsNVE) = sgMax_dummy(cellsNVE).*H(cellsNVE)./snr;
%        
%        if ~isempty(cellsNVEHorz) % NVE columns with horizontal flux from neighbors
%            % ---
%            % Determine ratio r of new CO2 to be distributed to top/bottom
%            r = value(opt.r);
%            % --- 
%             Sn = value(sG(cellsNVEHorz)); % current saturation
%             Snc = value(sgNVE(cellsNVEHorz)); % saturation (Sn < snr) at timestep when first nonzero flux from horizontal VE neighbors occured            
%             Sn_diff = max(Sn - Snc, 0);
%             %h_global(cellsNVEHorz) = (Sn - Snc).*H(cellsNVEHorz)./(1-swr); % Sn - Snc: additional saturation after NVE->VE transformation
%             h_global(cellsNVEHorz) = r.*Sn_diff.*H(cellsNVEHorz)./(1-swr);
%             %h_max_global(cellsNVEHorz) = Snc.*H(cellsNVEHorz)./snr + h_global(cellsNVEHorz);
%             h_max_global(cellsNVEHorz) = (1-r).*Sn_diff.*H(cellsNVEHorz)./snr + Snc.*H(cellsNVEHorz)./snr + value(h_global(cellsNVEHorz));
%        end   
%     end
%     
%     if ~isempty(cellsNVEMob)
%         % Distribute Sn as mobile CO2 at top and residual CO2 at bottom
%         pv = poreVolume(G, model.rock);
%         pv_mob = pv(cellsNVEMob);
%         h_global(cellsNVEMob) = (pv_mob.*sG(cellsNVEMob) - vG_B(cellsNVEMob))./(pv_mob.*(1-swr)).*H(cellsNVEMob);
%         h_max_global(cellsNVEMob) = vG_B(cellsNVEMob)./(pv_mob.*snr).*H(cellsNVEMob) + h_global(cellsNVEMob);
%     end
    % --------------------------   
    h_T_global = h_max_global.*(snr > 0) + h_global.*(snr == 0);
    h_B_global = zeros(size(sgMax));
    
    veTransition = op.connections.veToFineVertical | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;    
    veAll = op.connections.veInternalConn | op.connections.veTransitionHorizontalConn;
    
    n = op.N; 
    cn = op.N(veTransition, :);
    c_vic = op.N(veAll, :);
    veTrans = ismember(n, cn, 'rows');
    
    vG_bottom = vG_B(veTrans); % take abs val since fluxes are negative from bottom and up
    vG_bottom = abs(vG_bottom(veB));
    vG_bottom_smax = vG_B_smax(veTrans); % global indexing
    vG_bottom_smax = abs(vG_bottom_smax(veB)); % local indexing
    
    S = sG(cB); % VE bottom cells cB used directly as index because sG is ordered by cell numbers by default (1:G.cells.num)
    sMax = sgMax(cB).*(snr > 0) + S.*(snr == 0);
    
    S_B = 1./pv(cB) .* vG_bottom; % CO2 sat originating from bottom flux
    S_B_smax = 1./pv(cB) .* vG_bottom_smax; % earlier CO2 sat from bottom flux at time where smax was reached --> gives correct max saturation for top co2 part
    S_max_T = sMax - S_B_smax.*(snr > 0); % max saturation reached for top part (excluding part originating from bottom flux)
    
    h_T_global(cB) = max(H(cB).*S_max_T./(1-swr), 0); % max to avoid negative discrepancies
    h_T = h_T_global(cB); 
    % For h_B we use current saturation at bottom S_B, since
    % this is independent of smax
    h_B_global(cB) = max(h_T, H(cB).*(1-(S_B./(snr+eps)).*(snr>0))); % h_B will become h_T before h_T becomes negative (S_B > sMax), then h_max_global = H and col is in VE
    h_B = h_B_global(cB);              
    
    Snr_tot = (snr./H(cB)).*(h_T + (H(cB) - h_B)); % total residual saturation in column (NB: includes residual part in mobile zone!)
    S_mob = S - Snr_tot; % mobile saturation (at top of column)
       
    h_global(cB) = H(cB).*(S_mob./(1-swr-snr)); % subtract snr since this was included in Snr_tot
    h_max_global(cB) = H(cB) - h_B + h_T;
        
    % --------------------------
    
    p_entry = f.pe_rest;
    if any(vG_bottom)
        test = 0;
    end 
    
    %[pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h, H, rhoW, rhoG, muW, muG, isFine, isFine);
    [pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h_global, h_max_global, H, rhoW, rhoG, muW, muG, ...
                                                isFine, all_coarse_cells, h_global, h_max_global, {}, {}, ...
                                                cellsNVE, 've_model', opt.ve_model, 'p_entry', p_entry);
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
         upcw(transition_ve), upcg(transition_ve)] = computeTransitionFluxVE(model, pW, h_global, h_max_global, h_T_global, h_B_global, rhoW, rhoG, muW, muG, ...
                                                                              transition_ve, true, cellsNVE, cellsNVEMob, 'p_entry', p_entry);
    end
    % Treat transition between ve zones in vertical direction (for diffuse
    % leakage)
    transition_vertical = model.operators.connections.veToFineConn | ...
                          model.operators.connections.veTransitionVerticalConn & op.T > 0;
  
    if any(transition_vertical)        
        [vW(transition_vertical), vG(transition_vertical),...
         upcw(transition_vertical), upcg(transition_vertical)] = computeTransitionFluxVE(model, pW, h_global, h_max_global, h_T_global, h_B_global, rhoW, rhoG, muW, muG, ...
                                                                                            transition_vertical, false, cellsNVE, cellsNVEMob, 'p_entry', p_entry);
    end

end

function [pW, pW_f, sG, h, h_max, H, rhow, rhog, muw, mug, isFine] = getTransitionValuesVE_coarse(model, pW, h_global, h_max_global, h_T_global, h_B_global, ...
                                                                                                    index, subs, rhoW, rhoG, muW, muG, cellsNVE, cellsNVEMob, varargin)    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0);
    opt = merge_options(opt, varargin{:});
    
    c = model.operators.N(subs, index);
    cNVE = ismember(c, cellsNVE);
    cNVEMob = ismember(c, cellsNVEMob);
    
    isFine = model.G.cells.discretization(c) == 1;

    t = model.G.cells.topDepth(c);
    T = model.operators.connections.faceTopDepth(subs, index);
    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);

    %sG = getGasSaturationFromHeight(T, t, B, b, h_global(c));   
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
    %[a_M, a_R, sG] = getGasSatFromHeight(T, t, B, b, h_global(c), h_max_global(c), swr, snr, cNVE); 
    %[a_M, a_R, sG] = getGasSatFromHeightMob(T, t, B, b, h_global(c), h_max_global(c), swr, snr, cNVE, cNVEMob);
    [a_M, a_R, sG] = getGasSatFromHeightMob2(T, t, B, b, h_global(c), h_T_global(c), h_B_global(c), swr, snr);
    
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
%     pW = getPwFromHeight(B, t, b, h_global(c), h_max_global(c), pW(c), g, rhow, rhog, swr, snr, ...
%                             'cNVE', cNVE, 'p_entry', opt.p_entry, 'res_type', []);
    pW = getPwFromHeight(B, t, b, h_global(c), h_T_global(c), pW(c), g, rhow, rhog, swr, snr, ...
                                'cNVE', [], 'p_entry', opt.p_entry, 'res_type', []);
end

function [pW, pW_f, sG, h, h_max, H, rhow, rhog, muw, mug, isFine] = getTransitionValuesVE_fine(model, pW, h_global, h_max_global, h_T_global, h_B_global, index, subs, rhoW, rhoG, muW, muG, cellsNVE, cellsNVEMob, varargin)    
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0);
    opt = merge_options(opt, varargin{:});
    
    c = model.operators.N(subs, index);
    cNVE = ismember(c, cellsNVE); 
    cNVEMob = ismember(c, cellsNVEMob);
    t = model.G.cells.topDepth(c);
    
    isFine = model.G.cells.discretization(c) == 1;

    T = model.operators.connections.faceTopDepth(subs, index);

    b = model.G.cells.bottomDepth(c);
    B = model.operators.connections.faceBottomDepth(subs, index);
    H = model.operators.connections.faceHeight(subs, index);
    
    %sG = getGasSaturationFromHeight(T, t, B, b, h_global(c));
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);    
       
    %[a_M, a_R, sG] = getGasSatFromHeight(T, t, B, b, h_global(c), h_max_global(c), swr, snr, cNVE);   
    %[a_M, a_R, sG] = getGasSatFromHeightMob(T, t, B, b, h_global(c), h_max_global(c), swr, snr, cNVE, cNVEMob);    
    [a_M, a_R, sG] = getGasSatFromHeightMob2(T, t, B, b, h_global(c), h_T_global(c), h_B_global(c), swr, snr);   
    
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
%     pW_c = getPwFromHeight(C, t, b, h_global(c), h_max_global(c), pW_f, g, rhow, rhog, swr, snr, ...
%                             'cNVE', cNVE, 'p_entry', opt.p_entry, 'res_type', []);
    pW_c = getPwFromHeight(C, t, b, h_global(c), h_T_global(c), pW_f, g, rhow, rhog, swr, snr, ...
                            'cNVE', [], 'p_entry', opt.p_entry, 'res_type', []);
    pW = isFine.*pW_f + ~isFine.*pW_c;
end


function [vW, vG, upcw, upcg] = computeTransitionFluxVE(model, pW, h, h_max, h_T, h_B, rhoW, rhoG, muW, muG, vtc, treatAsCoarse, cellsNVE, cellsNVEMob, varargin)
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

        [pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesVE_coarse(model, pW, h, h_max, h_T, h_B, 1, vtc, rhoW, rhoG, muW, muG, cellsNVE, cellsNVEMob, 'p_entry', opt.p_entry);
        [pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesVE_coarse(model, pW, h, h_max, h_T, h_B, 2, vtc, rhoW, rhoG, muW, muG, cellsNVE, cellsNVEMob, 'p_entry', opt.p_entry);
      
        isFine_l = false;
        isFine_r = false;
        gdz_g = g*grad(t1, t2);
        gdz_w = g*grad(b1, b2);
    else
        % Treat coarse cells as fine cells
        % Original VE converted to fine cells => change h_max:
        
        [pW_l, pW_bl, sG_l, h_l, h_max_l, H_l, rhoW_l, rhoG_l, muW_l, muG_l] = getTransitionValuesVE_fine(model, pW, h, h_max, h_T, h_B, 1, vtc, rhoW, rhoG, muW, muG, cellsNVE, cellsNVEMob, 'p_entry', opt.p_entry);       
        [pW_r, pW_br, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r] = getTransitionValuesVE_fine(model, pW, h, h_max, h_T, h_B, 2, vtc, rhoW, rhoG, muW, muG, cellsNVE, cellsNVEMob, 'p_entry', opt.p_entry);

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
                                                         isFine_l, n1, h, h_max, vtc, 1, cellsNVE, ...
                                                         've_model', opt.ve_model, 'p_entry', opt.p_entry);
    [pW_r, pG_r, mobW_r, mobG_r] = evaluatePropertiesVE(model, pW_r, sG_r, h_r, h_max_r, H_r, rhoW_r, rhoG_r, muW_r, muG_r, ...
                                                         isFine_r, n2, h, h_max, vtc, 2, cellsNVE, ...
                                                         've_model', opt.ve_model, 'p_entry', opt.p_entry);    


    dpG   = grad(pG_l, pG_r) - rhoGf .* gdz_g;
    upcg  = (value(dpG)<=0);
    
    vG = -upstr(upcg, mobG_l, mobG_r) .* trans .* dpG;

    dpW   = grad(pW_l, pW_r) - rhoWf .* gdz_w;
    upcw  = (value(dpW)<=0);
    vW = -upstr(upcw, mobW_l, mobW_r) .* trans .* dpW;
end


function [pW, pG, mobW, mobG] = evaluatePropertiesVE(model, pW, sG, h, h_max, H, rhoW, rhoG, muW, muG, ...
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
    
end

function [pW, pG, mobW, mobG] = evaluatePropertiesVE2(model, pW, sG, h, h_T, h_B, H, rhoW, rhoG, muW, muG, ...
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

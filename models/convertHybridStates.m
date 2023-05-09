function states = convertHybridStates(model, model_fine, states_c, varargin)
% Converts hybrid states to states representative for the fine scale.
% Based on analytical reconstruction operators for pressure and saturation.

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

    ns = numel(states_c);
    states = cell(ns, 1);    
    states0_sGmax = states_c{1}.sGmax;   
    sgMax = states0_sGmax(model.G.partition);
    
    cellsBH = {states_c{1}.cBottom, states_c{1}.cHorz, states_c{1}.cBottomHorz};
    
    for i = 1:ns
        if i > 1
            states0_c = states_c{i-1};
        else
            states0_c = states_c{i};
        end
        [states{i}, sgMax] = convertState(model, model_fine, states_c{i}, states0_c, sgMax, cellsBH, varargin{:});        
        % what about converting flux?
    end
end

function [state, sgMax] = convertState(model, model_fine, state_c, state0_c, sgMax_fine, ...
                                                               cellsBH, varargin)
    % Convert a hybrid state to a fine state.
    % INPUTS:
    %   model: hybrid model to convert state for
    %   model_fune: full-dimensional model
    %   state_c: current hybrid state
    %   state0_c: previous hybrid state
    %   sgMax_fine: max saturation reached in fine-partitioned hybrid cells
    %   cellsBH: cell array of RVE cells in hybrid grid
    % RETURNS:
    %   state: reconstructed fine state
    %   sgMax: max gas saturation reached up to current state
    
    opt = struct('convert_flux', false, 'schedule', []);
    opt = merge_options(opt, varargin{:});
    
    f = model.fluid;
    ff = model_fine.fluid;
    op = model.operators;
    opf = model_fine.operators;

    CG = model.G;
    p = CG.partition;
    
    t = CG.cells.topDepth(p);
    T = CG.parent.cells.topDepth;
    
    b = CG.cells.bottomDepth(p);
    B = CG.parent.cells.bottomDepth;
    
    H = CG.cells.height(p); 
   
    % --- RESIDUAL SATURATION ---         
    swr = model.fluid.krPts.w(1);
    snr = model.fluid.krPts.g(1);
   
    %h_max = @(sgMax, H) H.*(sgMax./(1-swr));
    h = state_c.h;
    h_T = state_c.h_T;
    h_B = state_c.h_B;        
    
    cHorz = cellsBH{2}; % RVE cells with only horizontal fluxes from semi-perm layers
    hHi = state_c.hHi; % top depth of plume originaing from horizontal flux of semi-permeable layer
    Hi = state_c.Hi; % height from top of parent cell to semi-perm layer
    
    cBottomHorz = cellsBH{3}; % RVE cells with both horizontal and diffuse leakage up from semi-perm layers
    hBHi = state_c.hBHi; % top depth of plume originating from horizontal OR upward flux from semi-perm layer
    BHi = state_c.BHi;
      
    
    [a_M, a_R, sg] = getGasSatFromHeightFine(model, T, t, B, b, h(p), h_T(p), h_B(p), swr, snr, ...
                                             cHorz, cBottomHorz, hHi, hBHi, Hi, BHi);
    % ---------------------------------------------------------------            
    
    state = state_c;
    g = norm(model.gravity);
    if isfield(state_c, 'rho')
        rhow = state_c.rho(p, 1);
        rhog = state_c.rho(p, 2);
    else
        if isprop(model, 'EOSModel')
            rhow = model.fluid.rhoOS;
        else
            rhow = model.fluid.rhoWS;
        end
        rhog = model.fluid.rhoGS;
        rhog = repmat(rhog, model.G.cells.num, 1);
        rhow = repmat(rhow, model.G.cells.num, 1);
    end
    isFine = CG.cells.discretization == 1;
    
    % shift pressure from center to bottom
    state_c.pressure(~isFine) = state_c.pressure(~isFine) + g.*rhow(~isFine).*CG.cells.height(~isFine)/2;    
    
    cz = (T + B)/2;
    pressure = state_c.pressure(p); % fine-sclae pressures initially defined at bottom   
    rhow = rhow(p);
    rhog = rhog(p); 
    pW = getPwFromHeight(cz, t, b, h(p), h_T(p), pressure, g, rhow, rhog, swr, snr);
    % use original pressure from internal fine cells
    pW(isFine(p)) = pressure(isFine(p));
    
    state.pressure = pW;
    state.s = [1-sg, sg]; 
    state.sGmax = max(sg, sgMax_fine); % max of new sG and current max sG
    sgMax = state.sGmax;
    
    
    if opt.convert_flux
        % --- Probably not complete yet! ---

        vW = zeros(CG.parent.faces.num, 1);
        vG = zeros(CG.parent.faces.num, 1);
        % apply relperm function on reconstructed saturations
        krW = f.krW(1-sg); 
        krG = f.krG(sg);
        mobW = krW./f.muW(p);
        mobG = krG./f.muG(p);

        z = (T + B)/2; % midpoint of each cell in parent grid
        rhoWf = opf.faceAvg(rhow);
        dpW = opf.Grad(pW) - rhoWf.*g.*opf.Grad(z);
 
        sealingCells = find(ismember(p, CG.sealingCells));
        all_cells = (1:CG.parent.cells.num)';
        n_sealing = ismember(all_cells, sealingCells);
        pcap = ff.pcWG(sg);
        pG = pW + pcap; % all cells treated as fine
        rhoGf = opf.faceAvg(rhog);
        dpG = opf.Grad(pG) - rhoGf.*g.*opf.Grad(z);

        % parent faces
        faces = (1:CG.parent.faces.num)';
        bf = boundaryFaces(CG.parent); % boundary faces of parent gris
        ifaces = faces(~ismember(faces, bf)); % interior faces of parent grid
        open_bf = opt.schedule.control(1).bc.face;
        closed_bf = bf(~ismember(bf, open_bf));
        
        % Reconstruct trans over fine faces
        trans_p = model_fine.operators.T_all;    
        trans_int = trans_p(ifaces); % coarse trans uniformly distributed to associated fine cells               

        upcw  = (value(dpW)<=0);       
        upcg  = (value(dpG)<=0);
        mupW = opf.faceUpstr(upcw, mobW); % water mobility upwinding
        mupG = opf.faceUpstr(upcg, mobG);  % co2 mobility upwinding
               
        % map from faces to cells   
        obf2c = CG.parent.faces.neighbors(open_bf, :);
        obf2c = obf2c(obf2c ~= 0);
        % assign fluxes for interior faces
        vW(ifaces) = -mupW .* trans_int .* dpW;
        vG(ifaces) = -mupG .* trans_int .* dpG; 
        
        % Reconstruct for exterior faces       
        trans_ext = trans_p(open_bf); % transmissibility for boundary faces (uniform values by partitioning of fine grid)
        
        open_pW = opt.schedule.control(1).bc.value;             
        dpw = open_pW - pW(obf2c);
        open_pG = open_pW + pcap(obf2c); % add capillary pressure at open boundary
        dpg = open_pG - pG(obf2c);
        
        z_cell = (CG.parent.cells.bottomDepth + CG.parent.cells.topDepth)/2;        
        z_face = CG.parent.faces.centroids(:,3);    
        %gdz(open_bf) = g.*(z_face(open_bf) - z_cell(obf2c));    
        gdz = g.*(z_face(open_bf) - z_cell(obf2c));
        %dpG   = op.Grad(pG) - rhoGf .* gdz_g;
        dpG = dpg - rhog(obf2c).*gdz;
        dpW = dpw - rhow(obf2c).*gdz;
        
        % choose sign from non-zero neighbor cell at boundary interface
        sgn = 1 - 2*(CG.parent.faces.neighbors(open_bf, 2) == 0); % -1 (if first neighbor is 0) or 1 (if second neighbor is 0)
        mupG = sgn.*mobG(obf2c);
        mupW = sgn.*mobW(obf2c); %.*(value(dpW)<=0);                      
        
        vG(open_bf) = -mupG.*trans_ext.*dpG; % bf2c is always upwind cell for half-faces
        vG(closed_bf) = 0;
        vW(open_bf) = -mupW.*trans_ext.*dpW;
        vW(closed_bf) = 0;
 
        state.flux = zeros(CG.parent.faces.num, 2);
        state.flux(:,1) = vW;
        state.flux(:,2) = vG;
        
        % --- OR USE: ---
        %getBoundaryConditionFluxesAD by inputting reconstructed pressure
        %and sat
        % ---------------
    end

    if isprop(model, 'EOSModel')
        state.T = state.T(p);
        state.L = state.L(p);
        state.K = state.K(p, :);
        state.components = state_c.x(CG.partition, :).*state.s(:, 1) + state_c.y(CG.partition, :).*state.s(:, 2);
        pureLiquid = state.s(:, 1) == 1;
        pureVapor = state.s(:, 2) == 1;
        
        state.x = ~pureLiquid.*state_c.x(CG.partition, :) + pureLiquid.*state.components;
        state.y = ~pureVapor.*state_c.y(CG.partition, :) + pureVapor.*state.components;
    end
end


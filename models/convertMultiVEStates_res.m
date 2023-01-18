function states = convertMultiVEStates_res(model_c, model_f, states_c, states_f, varargin)
%Undocumented Utility Function

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
    sgMax = states0_sGmax(model_c.G.partition);
    snMaxVE_bottom = 0;
    vG_any = 0; % only CO2 flux at initial state is at well
    
    for i = 1:ns        
        [states{i}, sgMax, snMaxVE_bottom, vG_any] = convertState(model_c, model_f, states_c{i}, states_f{i}, sgMax, snMaxVE_bottom, vG_any, i, varargin{:});        
        %states0_sGmax = states{i}.sGmax; % update historically max fine-scale sat
        % what about converting flux?
    end
end

function [state, sgMax, snMaxVE_bottom, vG_any] = convertState(model_c, model_f, state_c, state_f, sgMax_fine, snMaxVE_bottom, vG_any, i, varargin)
    opt = struct('convert_flux', false, 'schedule', []);
    opt = merge_options(opt, varargin{:});
    
    f = model_c.fluid;
    op = model_c.operators;

    CG = model_c.G;
    p = CG.partition;
    
    all_cells = ismember(p, (1:CG.cells.num)');    
    
    t = CG.cells.topDepth(p);
    T = CG.parent.cells.topDepth;   
    b = CG.cells.bottomDepth(p);
    B = CG.parent.cells.bottomDepth;    
    H = CG.cells.height(p); 
   
    % --- RESIDUAL SATURATION ---         
    swr = model_c.fluid.krPts.w(1);
    snr = model_c.fluid.krPts.g(1);

    sgMax = state_c.sGmax;
    sgMax_c = sgMax(p); % max saturation from hybrid - not updated since coarse states already solved for
    sG = state_c.s(:,2);
    sG_c = sG(p);
      
    vG = state_c.vGsum; % get max absolute co2 flux for each coarse connection
       
    h_max_func = @(sgMax, H) H.*(sgMax./(1-swr));   
    h_func = @(sg, sgMax, H) H.*((sg.*(1-swr) - sgMax.*snr) ./ ((1-swr).*(1-swr-snr)));  
    
    [sg, h, h_max, cellsNVE] = height2SatConvert(model_c, h_func, h_max_func, sG, sgMax, vG, i);
    % Returned h and h_max are for parent grid, same value copied for all
    % fine cells of a VE column.
    cellsNVE = ismember(p, cellsNVE);
    % ---------------------------------------------------------------            
    
    state = state_c;
    g = norm(model_c.gravity);
    if isfield(state_c, 'rho')
        rhow = state_c.rho(p, 1);
        rhog = state_c.rho(p, 2);
    else
        if isprop(model_c, 'EOSModel')
            rhow = model_c.fluid.rhoOS;
        else
            rhow = model_c.fluid.rhoWS;
        end
        rhog = model_c.fluid.rhoGS;
        rhog = repmat(rhog, model_c.G.cells.num, 1);
        rhow = repmat(rhow, model_c.G.cells.num, 1);
    end
    isFine = CG.cells.discretization == 1;
    
    state.s = [1-sg, sg]; 
    % choose max among current calculated fine scale sat and fine scale sat
    % from all previous states
    state.sGmax = max(sg, sgMax_fine); % max of new sG and current max sG
    sgMax = state.sGmax;
    
    % shift pressure from center to bottom
    state_c.pressure(~isFine) = state_c.pressure(~isFine) + g.*rhow(~isFine).*CG.cells.height(~isFine)/2;    
    % shift pressure from center to top of each VE column
    %state_c.pressure(~isFine) = state_c.pressure(~isFine) - g.*rhow(~isFine).*CG.cells.height(~isFine)/2;
    
    cz = (T + B)/2;
    pressure = state_c.pressure(p); % fine-sclae pressures initially defined at top
    % WHY USE THIS FUNCTION TO RECONSTRUCT PRESSURE IN INTERNAL VE COLUMNS? DOESN'T DUPUIT ASSUMPTION HOLD HERE?  
    % Dupuit not valid only for veToFine and veVerticalConn -> significant
    % vertical flow here (captured by fine cells)
    rhow = rhow(p);
    rhog = rhog(p);
    %p_c = getPwFromHeight(cz, t, b, h(sG_c, sgMax_c, H), h_max(sgMax_c, H), pressure, g, rhow, rhog, swr, snr);    
    %p_c = getPwFromHeight(cz, t, b, h(sg, sgMax, H), h_max(sgMax, H), pressure, g, rhow, rhog, swr, snr, cellsNVE, all_cells);
    cNVE = ismember(all_cells, cellsNVE);
    p_c = getPwFromHeight(cz, t, b, h, h_max, pressure, g, rhow, rhog, swr, snr, 'cNVE', cNVE);
    % use original pressure from internal fine cells
    p_c(isFine(p)) = pressure(isFine(p));
    % use Dupuit approx in internal ve conn (WHAT ABOUT HORIZONTAL VE
    % CONN -> can't assume Dupuit hold here?)
    vIc_cells = zeros(CG.cells.num, 1); % veInternalConn cells
    vIc_n1 = op.N(op.connections.veInternalConn, 1);    
    vIc_n2 = op.N(op.connections.veInternalConn, 2);
    vIc_cells(vIc_n1 | vIc_n2) = 1;
    dupuit = logical(vIc_cells(p));   
    p_c(dupuit) = pressure(dupuit) - rhow(dupuit).*g.*cz(dupuit); % pressure from bottom back to center
    
    state.pressure = p_c;    
    
    % -------- Reconstruct fluxes ------------
    if opt.convert_flux
        vW = zeros(CG.parent.faces.num, 1);
        vG = zeros(CG.parent.faces.num, 1);
        % apply relperm function on reconstructed saturations
        krW = f.krW(1-sG); 
        krG = f.krG(sG);
        mobW = krW./f.muW;
        mobG = krG./f.muG;

        z = (T + B)/2; % midpoint of each cell in parent grid
        dpW = op.Grad(p_c) - rhow.*g.*op.Grad(z);

        sealingCells = CG.sealingCells(p); % sealing cells for parent grid
        all_cells = (1:CG.parent.cells.num)';
        n_sealing = ismember(all_cells, sealingCells);
        pcap = f.pcWG(sG, n_sealing);
        pG = p_c + pcap; % all cells treated as fine
        dpG = op.Grad(pG) - rhog.*g.*op.Grad(z);

        % Reconstruct trans over fine faces
        trans_c = model_c.operators.T_all;
        trans_p = trans_c(CG.faces.fconn);
        trans_int = trans_p(bf); % coarse trans uniformly distributed to associated fine cells               

        upcw  = (value(dpW)<=0);       
        upcg  = (value(dpG)<=0);
        
        % parent faces
        faces = (1:CG.parent.faces.num)';
        bf = boundaryFaces(CG.parent); % boundary faces of parent gris
        ifaces = faces(~ismember(faces, bf)); % interior faces of parent grid
        open_bf = opt.schedule.control(1).face;
        closed_bf = bf(~ismember(bf, open_bf));
        % map from faces to cells
        if2c = CG.parent.faces.neighbors(ifaces, :);
        bf2c = CG.parent.faces.neighbors(bf, :);
        bf2c = bf2c(bf2c ~= 0); % choose cell at boundary but inside grid    
        obf2c = CG.parent.faces.neighbors(open_bf, :);
        obf2c = obf2c(obf2c ~= 0);
        % assign fluxes for interior faces
        vW(ifaces) = -op.faceUpstr(upcw, mobW(if2c)) .* trans_int .* dpW(ifaces);
        vG(ifaces) = -op.faceUpstr(upcg, mobG(if2c)) .* trans_int .* dpG(ifaces); 
        
        % Reconstruct for exterior faces
        trans_c = model_c.operators.T_all;
        trans_p = trans_c(CG.faces.fconn);
        trans_ext = trans_p(bf); % transmissibility for boundary faces (uniform values by partitioning of fine grid)
        
        open_pW = opt.schedule.control(1).value;       
        dpw = op.Grad(pW);
        dpw(open_bf) = open_pW - pW(obf2c); % boundary face water pressure - boundary cell water pressure
        dpw(closed_bf) = 0;
        open_pG = open_pW + pcap(obf2c); % add capillary pressure at open boundary
        dpg = op.Grad(pG);
        dpg(open_bf) = open_pG - pG(obf2c); % boundary face gas pressure - boundary cell gas pressure
        dpg(closed_bf) = 0; 
        
        z_cell = (CG.parent.cells.bottomDepth + CG.parent.cells.topDepth)/2;
        gdz = g.*z_cell;
        z_face = CG.parent.faces.centroids(:,3);    
        %gdz(open_bf) = g.*(z_face(open_bf) - z_cell(obf2c));    
        gdz(bf) = g.*(z_face(bf) - z_cell(bf2c));
        %dpG   = op.Grad(pG) - rhoGf .* gdz_g;
        dpG = dpg(bf) - rhoGf.*gdz(bf);
        dpW = dpw(bf) - rhoWf.*gdz(bf);
        
        vG(bf) = mobG(bf2c).*trans_ext.*dpG; % bf2c is always upwind cell for half-faces
        vG(closed_bf) = 0;
        vW(bf) = mobW(bf2c).*trans_ext.*dpW;
        vW(closed_bf) = 0;
                
        %cts = rldecode((1:CG.faces.num)', diff(CG.faces.connPos));
        %Tc = accumarray(cts, model_f.operators.T_all(CG.faces.fconn));
        state.flux = zeros(CG.parent.faces.num, 2);
        state.flux(:,1) = vW;
        state.flux(:,2) = vG;
    end
    % ---------------------------------------
    
    if isprop(model_c, 'EOSModel')
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

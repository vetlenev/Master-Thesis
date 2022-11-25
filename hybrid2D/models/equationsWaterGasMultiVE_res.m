function [problem, state] = equationsWaterGasMultiVE_res(model, state0, state, dt, drivingForces, varargin)
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

    opt = struct('Verbose'     , mrstVerbose , ...
                 'reverseMode' , false       , ...
                 'resOnly'     , false       , ...
                 'iteration'   , -1          , ...
                 'stepOptions' , []); % compatibility only
    opt = merge_options(opt, varargin{:});
    
    assert(isempty(drivingForces.src));  % unsupported
    op = model.operators;
    G = model.G;
    p = G.partition;
    f = model.fluid;
    
    % Extract current and previous values of all variables to solve for
    [pW, sG, wellSol] = model.getProps(state, 'pressure', 'sg', 'wellsol');   
    [p0, sG0, wellSol0] = model.getProps(state0, 'pressure', 'sg', 'wellsol');
    
    [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

    % ------------------ Initialization of independent variables ------------------
    backend = model.AutoDiffBackend;
    if ~opt.resOnly
        if ~opt.reverseMode
            [pW, sG, wellVars{:}] = backend.initVariablesAD(pW, sG, wellVars{:});
        else
            wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
            [p0, sG0, wellVars0{:}] = backend.initVariablesAD(p0, sw0, wellVars0{:}); %#ok
        end
    end
    
    % Set up AD state
    sat = {1-sG, sG};
    sat0 = {1-sG0, sG0};
    state = model.setProps(state, {'s', 'pressure'}, {sat, pW});
    state0 = model.setProps(state0, {'s', 'pressure'}, {sat0, p0});
    % Set up properties
    state = model.initStateFunctionContainers(state);
    
    % Get properties etc
    [b, pv, rho, mu, p_phase] = model.getProps(state, 'ShrinkageFactors', 'PoreVolume', 'Density', 'Viscosity', 'PhasePressures');
    [b0, pv0] = model.getProps(state0, 'ShrinkageFactors', 'PoreVolume');
    [bW, bG] = deal(b{:});
    [muW, muG] = deal(mu{:});
    [bW0, bG0] = deal(b0{:});
    
    trans = model.getProps(state, 'Transmissibility');
    sgMax = model.getProps(state, 'sGmax');    

    g = norm(model.gravity);
    pW_center = pW;
    
    isVE = G.cells.discretization > 1;  
    % ----- Find relevant VE bottom cells
    veTransition = op.connections.veToFineVertical | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;        
    veAll = op.connections.veInternalConn | op.connections.veTransitionHorizontalConn;
    veHorz = op.connections.veTransitionHorizontalConn;
      
    c_bottom = op.N(veTransition, :);
    c_horz = op.N(veHorz, :); 
    c_vic = op.N(veAll, :);
    cB = {}; veB = {};
    cH = {}; veH = {};           
   
    for idx=1:2
        % treat bottom fluxes
        c = c_bottom(:,idx);
        isVE_c = isVE(c);
        if any(isVE_c)
            t = G.cells.topDepth(c);
            T = op.connections.faceTopDepth(veTransition, idx);
            b = G.cells.bottomDepth(c);
            B = op.connections.faceBottomDepth(veTransition, idx);
            Hb = G.cells.height(c); % to not overwrite global H
            
            veb = B == b & T ~= t; 
            cb = c(veb); % select correct bottom cells 
            %cb = c_vic(ismember(c_vic, cb));            
            if ~isempty(veb)
                veB = cat(1, veB, find(veb));
                cB = cat(1, cB, cb);
            end
        end
         
        % treat horizontal fluxes (all cells connected to veHorizontal are
        % VE)
        c = c_horz(:,idx); % horizontal cells in coarse grid
        %cp = pc_horz(:,idx); % horizontal cells in parent grid
        
        t = G.cells.topDepth(c);
        tp = G.cells.topDepth(p);
        %tp = tp(cp);
        T = op.connections.faceTopDepth(veHorz, idx);     
        b = G.cells.bottomDepth(c);
        B = op.connections.faceBottomDepth(veHorz, idx);
        Hb = G.cells.height(c); % to not overwrite global H

        veH{idx} = find(T ~= t); % all virtual ve cells except top one (this does not give bottom saturation region)
        cH{idx} = c(veH{idx}); % select correct ve virtual cells
        %cV = c_vic(ismember(c_vic, cv));        
    end
    %cV = unique(cell2mat(cV)); 
    [cB, cB_idx, ~] = unique(cell2mat(cB));
    veB = cell2mat(veB);
    veB = veB(cB_idx); % <----- shuffle connections in same order as unique cells
    % ----- 
      
    % Modify the fluxes to account for VE transitions
    [rhoW, rhoG] = deal(rho{:});
    
    % For VE cells, water pressure (primary var) evaluated at bottom
    pW = pW_center + isVE.*rhoW.*g.*G.cells.height/2; % center -> bottom    
    %pW = pW_center - isVE.*rhoW.*g.*G.cells.height/2; %  center -> top       
    
    % --- Accounting for residual sat in coarse model ---       
    swr = f.krPts.w(1);   
    snr = f.krPts.g(1);
    
    if ~isfield(state0, 'vGsum')
        state0.vGsum = zeros(size(model.operators.N, 1), 1); % initially zero fluxes by default
        %state0.vGsum = zeros(numel(cB), 1);
    end
    if ~isfield(state0, 'sGnve')
       state0.sGnve = cell2mat(state0.s(:,2)); % initial saturation
    end
    if ~isfield(state0, 'vGsMax')
       state0.vGsMax = zeros(size(model.operators.N, 1), 1);
       state.vGsMax = state0.vGsMax;
       %state0.vGsMax = zeros(numel(cB), 1); 
    end
    if ~isfield(state0, 'h_T') || ~isfield(state0, 'h_B') || ~isfield(state0, 'h')
       state0.h_T = zeros(size(state0.sGnve)); 
       state0.h_B = zeros(size(state0.sGnve));
       state0.h = zeros(size(state0.sGnve));
    end
       
    % Get cells partly residual filled and fully residual filled *from below*
    [c_prf, c_frf] = getResidualFilledCells(model, sG, state0.vGsum); % CHANGED FROM sG to sgMax !!
    %[c_prf, c_mrf, c_frf] = getResidualFilledCellsMob(model, pv0, bG0, sG0, sgMax, state0.vGsum); % [partly residual filled, mobile residual filled, fully residual filled]    
    
    %sg = value(sG); % to avoid ADI/double warning when using sG for caluclations   
    [vW, vG, mobW, mobG, upcw, upcg] = computeHybridFluxesVEres(model, pW, sG, muW, muG, rhoW, rhoG, trans, sgMax, c_prf, c_frf);  
    %[vW, vG, mobW, mobG, upcw, upcg, ...
    %    h, h_T, h_B, hHi_state, hBHi_state] = computeHybridFluxesVEres_test(model, pW, sG, muW, muG, rhoW, rhoG, trans, sgMax, state0.vGsum, state0.vGsMax, cB, veB, cH, veH);
            
%     state.h = h;
%     state.h_T = h_T;
%     state.h_B = h_B;
%     
%     state.cHorz = hHi_state{1}; % VE horizontal transition cells fulfilling T ~= t
%     state.hHi = hHi_state{2}; % depth of (top of) bottom plume in each virtual VE cell
%     state.Hi = hHi_state{3}; % depth of bottom of associated virtual VE cells
%     
%     state.cBottomHorz = hBHi_state{1};
%     state.hBHi = hBHi_state{2};
%     state.BHi = hBHi_state{3};
    
    %state.vGsum = max(abs(vG), state0.vGsum);   
    % ---------------------------------------- -----------
    
    bWvW = op.faceUpstr(upcw, bW).*vW;
    bGvG = op.faceUpstr(upcg, bG).*vG;
     
    state.vGsum = state0.vGsum + bGvG.*dt; % abs-value: assume all fluxes through bottom interface is directed upwards
       
    % --- BOTTOM flux summed up to time step of first occurence of current sgMax
%     veBottom = ismember(op.N, c_bottom, 'rows');
%     vGsum_bottom = value(state.vGsum(veBottom));
%     vGsMax_bottom = state0.vGsMax(veBottom); % NB: important to choose from earlier state
% 
%     veB_global = find(veBottom);
%     veB_global = veB_global(veB); % global index connection (veB is local index connection for veToFine and veVertical transition connections)
%     state.vGsMax(veB_global) = vGsMax_bottom(veB).*(value(sG(cB)) < sgMax(cB)) + ... % choose stored summed bottom flux
%                                 vGsum_bottom(veB).*(value(sG(cB)) >= sgMax(cB)); % choose current summed bottom flux
%                             
    % --- HORIZONTAL fluxes summed up to time step of first occurence of current sgMax
%     veHorz = ismember(op.N, c_horz, 'rows');                        
%     vGsum_horz = value(state.vGsum(veBottom)); % Choose current summed flux (since sG exceeds sMax)
%     vGsMax_horz = state0.vGsMax(veBottom); % Choose from earlier state (point when sMax was reached)
% 
%     veH_global = find(veHorz);
%     veH = [veH{1}; veH{2}];
%     cH = [cH{1}; cH{2}];
%     veH_global = veH_global(veH); % global index connection (veB is local index connection for veToFine and veVertical transition connections)
%     state.vGsMax(veH_global) = vGsMax_horz(veH).*(value(sG(cH)) < sgMax(cH)) + ... % choose stored summed bottom flux
%                                 vGsum_horz(veH).*(value(sG(cH)) >= sgMax(cH)); % choose current summed bottom flux
%                             
    % ---------------
    
    if model.outputFluxes
        state = model.storeFluxes(state, vW, [], vG);
    end
    if model.extraStateOutput  
        state = model.storebfactors(state, bW, [], bG);
        state = model.storeMobilities(state, mobW, [], mobG);
        state = model.storeUpstreamIndices(state, upcw, [], upcg);
    end
      
    % --------------------------- Continuity equations ---------------------------
     
    eqs{1} = (1/dt) .* (pv .* bW .* (1-sG) - pv0 .* bW0 .* (1-sG0)) + op.Div(bWvW);
    eqs{2} = (1/dt) .* (pv .* bG .* sG     - pv0 .* bG0 .* sG0    ) + op.Div(bGvG);

    % ------------------------------ Well equations ------------------------------
    mob = {mobW, mobG};

    primaryVars = ['pressure' , 'sG', wellVarNames];
    types = {'cell'           , 'cell'};
    names = {'water'          , 'gas'};

    [eqs, state] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                     p_phase, sat, mob, rho, ...
                                                                     {}, {}, ...
                                                                     drivingForces);
    
    
    [eqs, names, types, state.wellSol] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, pW_center, mob, rho, {}, {}, dt, opt);

    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);    
    % END
end




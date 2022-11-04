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
    
    % Modify the fluxes to account for VE transitions
    [rhoW, rhoG] = deal(rho{:});
   
    % For VE cells, water pressure (primary var) evaluated at bottom
    pW = pW_center + isVE.*rhoW.*g.*G.cells.height/2; % center -> bottom    
    %pW = pW_center - isVE.*rhoW.*g.*G.cells.height/2; %  center -> top       
    
    % --- Accounting for residual sat in coarse model ---       
    swr = f.krPts.w(1);   
    snr = f.krPts.g(1);
    
    if ~isfield(state0, 'vGmax')
        state0.vGmax = zeros(size(model.operators.N, 1), 1); % initially zero fluxes by default
    end
    if ~isfield(state0, 'sGnve')
       state0.sGnve = cell2mat(state0.s(:,2)); 
    end
       
    % Get cells partly residual filled and fully residual filled *from below*
    [c_prf, c_frf] = getResidualFilledCells(model, sG, state0.vGmax); % CHANGED FROM sG to sgMax !!
    %[c_prf, c_frf, c_horz] = getResidualFilledCells_test(model, sG, state0.vGmax, sG0, state0.sGnve);   
        
    % NB: Elements of sgMax corresponding to partially residual filled cells
    % are not updated in the state-object because we manually modify these 
    % elements of sgMax but not sG. But saturation for all other cells 
    % (i.e. non-NVE) are updated in the state-object.
    state.sGnve = sG; % store current saturation
    
%     if ~isempty(c_horz) % retrieve saturation from time step when transitioning from NVE to VE %state0.sGnve == 0 && state.sGnve > 0
%         state.sGnve(c_horz) = state0.sGnve(c_horz);           
%     end

    %state = model.setProps(state, 'sGmax', sgMax);   
    
    [vW, vG, mobW, mobG, upcw, upcg] = computeHybridFluxesVEres(model, pW, sG, muW, muG, rhoW, rhoG, trans, sgMax, c_prf, c_frf);  
    %[vW, vG, mobW, mobG, upcw, upcg] = computeHybridFluxesVEres_test(model, pW, sG, muW, muG, rhoW, rhoG, trans, sgMax, state.sGnve, c_prf, c_frf, c_horz);
       
    state.vGmax = max(abs(vG), state0.vGmax);
    if ~isa(sG, 'ADI')
       stop = 0; 
    end
    % ---------------------------------------- -----------
    
    bWvW = op.faceUpstr(upcw, bW).*vW;
    bGvG = op.faceUpstr(upcg, bG).*vG;
    
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
end




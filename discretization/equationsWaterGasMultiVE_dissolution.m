function [problem, state] = equationsWaterGasMultiVE_dissolution(model, state0, state, dt, drivingForces, varargin)
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
                 'schedule', [], ...
                 'stepOptions' , []); % compatibility only
    opt = merge_options(opt, varargin{:});
    
    assert(isempty(drivingForces.src));  % unsupported
    op = model.operators;
    G = model.G;    
    p = G.partition;
    f = model.fluid;   
    
    % Extract current and previous values of all variables to solve for
    [pW, sG, sGmax, rs, wellSol] = model.getProps(state, 'pressure', 'sg', 'sGmax', 'rs', 'wellsol');   
    [p0, sG0, sGmax0, rs0, wellSol0] = model.getProps(state0, 'pressure', 'sg', 'sGmax', 'rs', 'wellsol');
    
    [wellVars, wellVarNames, wellMap] = model.FacilityModel.getAllPrimaryVariables(wellSol);

    % All cells with CO2-saturated brine
    rsSat  = f.rsSat(pW);
    s_tol  = 0; %(f.dis_rate == 0) * sqrt(eps); % small nonzero if instantaneous dissolution model
    isSat  = (sG > s_tol) | (rs > rsSat);

    % ------------------ Initialization of independent variables ------------------
    backend = model.AutoDiffBackend;
    if ~opt.resOnly
        if ~opt.reverseMode
            [pW, sG, sGmax, rs, wellVars{:}] = backend.initVariablesAD(pW, sG, sGmax, rs, wellVars{:});
        else
            wellVars0 = model.FacilityModel.getAllPrimaryVariables(wellSol0);
            [p0, sG0, sGmax0, rs0, wellVars0{:}] = backend.initVariablesAD(p0, sG0, sGmax0, rs0, wellVars0{:}); %#ok
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
    
    [pvMult, transMult, mobMult, pvMult0] = getMultipliers(f, pW, p0);

    g = norm(model.gravity);
    pW_center = pW;
    
    isVE = G.cells.discretization > 1;  
    % ----- Find relevant VE bottom cells
    veTransition = op.connections.veToFineVertical | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;        
    veHorz = op.connections.veTransitionHorizontalConn;
      
    c_bottom = op.N(veTransition, :);
    c_horz = op.N(veHorz, :);     
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
            
            veb = B == b & T ~= t; % query only semi-perm layers at bottom of VE col 
            cb = c(veb); % select correct bottom cells                        
            if ~isempty(veb)
                veB = cat(1, veB, find(veb));
                cB = cat(1, cB, cb);
            end
        end
         
        % treat horizontal fluxes (all cells connected to veHorizontal are
        % VE => no need to check any(isVE_c))
        c = c_horz(:,idx); % horizontal cells in coarse grid
        
        t = G.cells.topDepth(c);        
        T = op.connections.faceTopDepth(veHorz, idx);     
        
        veH{idx} = find(T ~= t); % all virtual ve cells except top one (this does not give bottom saturation region)
        cH{idx} = c(veH{idx}); % select correct ve virtual cells   
    end

    [cB, cB_idx, ~] = unique(cell2mat(cB));
    veB = cell2mat(veB);
    veB = veB(cB_idx); % <----- shuffle connections in same order as unique cells
    
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
    if ~isfield(state0, 'vGsum_sMax') % to store accumulated flux up to point where current sMax
       state0.vGsum_sMax = zeros(size(model.operators.N, 1), 1);
       state.vGsum_sMax = state0.vGsum_sMax; 
    end
    if ~isfield(state0, 'h_T') || ~isfield(state0, 'h_B') || ~isfield(state0, 'h')
       state0.h_T = zeros(size(state0.sGnve)); 
       state0.h_B = zeros(size(state0.sGnve));
       state0.h = zeros(size(state0.sGnve));
    end
           
    [vW, vG, mobW, mobG, upcw, upcg, ...
        h, h_max, h_T, h_B, cellsBH, hHi_state, hBHi_state] = computeHybridFluxesVE_test(model, pW, sG, muW, muG, rhoW, rhoG, trans, sGmax, state0.vGsum, state0.vGsum_sMax, cB, veB, cH, veH);
            
    state.h = h;
    state.h_max = h_max;
    state.h_T = h_T;
    state.h_B = h_B;
    % NB: Change to only append store cells ONCE (no need to store the same
    % cells for every state.
    if ~isfield(state0, 'cBottom') % only append (static) cells to one state to avoid redundant stores in later states
        state.cBottom = cellsBH{1}; state.fBottom = cellsBH{4}; % bottom cells; bottom connections
        state.cHorz = cellsBH{2}; state.fHorz = cellsBH{5}; % VE horizontal transition cells fulfilling T ~= t   
        state.cBottomHorz = cellsBH{3}; state.fBottomHorz = cellsBH{6};
    end
    
    state.hHi = hHi_state{1}; % depth of (top of) bottom plume (from horizontal flux) in each virtual VE cell
    state.Hi = hHi_state{2}; % depth of bottom of associated virtual VE cells
        
    state.hBHi = hBHi_state{1};
    state.BHi = hBHi_state{2};
    
    %state.vGsum = max(abs(vG), state0.vGsum);   
    % ---------------------------------------- -----------
    
    bWvW = op.faceUpstr(upcw, bW) .* vW;
    bGvG = op.faceUpstr(upcg, bG) .* vG;
    rsbWvW = op.faceUpstr(upcw, rs) .* bWvW;
     
    %state.vGsum = state0.vGsum + (bGvG + rsbWvW).*dt; % dissolved CO2 is
    %not part of gas phase, but we want accumulated flux of gas phase, not
    %CO2 components, so we should not include rsbWvW in sum.
    state.vGsum = state0.vGsum + bGvG.*dt;
       
    % --- BOTTOM flux summed up to time step of first occurence of current sGmax
    veBottom = ismember(op.N, c_bottom, 'rows');
    vGsum_bottom = value(state.vGsum(veBottom)); % accumulated flux from bottom, up to current time step
    vGsum_sMax_bottom = state0.vGsum_sMax(veBottom); % store accumulated flux from previous time step

    veB_global = find(veBottom);
    veB_global = veB_global(veB); % global index connection (veB is local index connection for veToFine and veVertical transition connections)
    state.vGsum_sMax(veB_global) = vGsum_sMax_bottom(veB).*(value(sG(cB)) < sGmax(cB)) + ... % choose stored summed bottom flux
                                vGsum_bottom(veB).*(value(sG(cB)) >= sGmax(cB)); % choose current summed bottom flux
    % state.vGsum_sMax == state0.vGsum_sMax if no new sMax reached this time step
                            
    % --- HORIZONTAL fluxes summed up to time step of first occurence of current sGmax
    veHorz = ismember(op.N, c_horz, 'rows');                        
    vGsum_horz = value(state.vGsum(veHorz));
    vGsum_sMax_horz = state0.vGsum_sMax(veHorz);

    veH_global = find(veHorz);
    veH = [veH{1}; veH{2}];
    cH = [cH{1}; cH{2}];
    veH_global = veH_global(veH); % global index connection (veB is local index connection for veToFine and veVertical transition connections)
    state.vGsum_sMax(veH_global) = vGsum_sMax_horz(veH).*(value(sG(cH)) < sGmax(cH)) + ... % choose stored summed bottom flux
                                vGsum_horz(veH).*(value(sG(cH)) >= sGmax(cH)); % choose current summed bottom flux
                             
    
    if model.outputFluxes
        state = model.storeFluxes(state, vW, [], vG);
    end
    if model.extraStateOutput  
        state = model.storebfactors(state, bW, [], bG);
        state = model.storeMobilities(state, mobW, [], mobG);
        state = model.storeUpstreamIndices(state, upcw, [], upcg);
    end
      
    sW = 1 - sG;
    sW0 = 1 - sG0;
    eqs = cell(1, 4);
    % --------------------------- Continuity equations ---------------------------   
    eqs{1} = (op.pv/dt) .* (pvMult .* bW .* sW - pvMult0 .* bW0 .* sW0) + op.Div(bWvW);
    eqs{2} = (op.pv/dt) .* (pvMult .* (bG.*sG + rs.*bW.*sW) ...
                        - pvMult0 .* (bG0.*sG0 + rs0.*bW0.*sW0)) ...
                     + op.Div(bGvG + rsbWvW);

    % ------------------------------ Well equations ------------------------------
    mob = {mobW, mobG};
    dissolved = {{[], []}, {rs, []}};

    primaryVars = ['pressure' , 'sG', 'rs', 'sGmax', wellVarNames];
    types = {'cell'           , 'cell', 'cell', 'cell'};
    names = {'water'          , 'gas', 'dissol', 'sGmax'};

    [eqs, state, src] = addBoundaryConditionsAndSources(model, eqs, names, types, state, ...
                                                                     p_phase, sat, mob, rho, ...
                                                                     dissolved, {}, ...
                                                                     drivingForces);       
    
    [eqs, names, types, state.wellSol, qWell] = model.insertWellEquations(eqs, names, types, wellSol0, wellSol, wellVars, wellMap, pW_center, mob, rho, dissolved, {}, dt, opt);
    wc = vertcat(drivingForces.W.cells);

    if isempty(wellMap.isRate)
        qW_pV = 0; % no phase volume in undefined wells
    else
        qW_pV = qWell.phaseVolume{1};
    end
    % Equations for dissolution    
    eqs(3:4) = compute_dissolution_equations(model, state, G, f, sG, sG0, sGmax, ...
                                            sGmax0, bG, sW, sW0, bW, bW0, pW_center, ...
                                            rhoW, mobW, rs, rs0, rsSat, isSat, ...
                                            pvMult, pvMult0, wc, qW_pV, op, ... % pvMult, pvMult0
                                            rsbWvW, src, dt, isVE);

    % Rescaling non-well equations
    for i = 1:4
       eqs{i} = eqs{i} * dt / year;
    end

    problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);    
    % END
end

function eqs = compute_dissolution_equations(model, state, G, f, sG, sG0, sGmax, ...
                                             sGmax0, bG, sW, sW0, bW, bW0, p, ...
                                             rhoW, mobW, rs, rs0, rsSat, isSat, ...
                                             pvMult, pvMult0, wc, cqsw, op, ...
                                             rsbWvW, src, dt, isVE)                 
      
      
%       eqs{1} = sG + 0*sG0;
%       eqs{1}(isSat) = rs(isSat) - rsSat(isSat); % set dissolved saturation rs to maximum possible (rsSat) for all cells with any CO2 (sG > 0) -> this CO2 will dissolve to the maximum possible instantaneously
%       eqs{2}        = sGmax - max(sGmax0, sG); % standard equation for sMax -> no dependency on time-derivative (as is the case for rate-driven equations) since changes happen instantaneously

      
      % VE cells: Rate-driven dissolution model       
      rate = f.dis_rate .* op.pv./G.cells.height; % rate per area multiplied
                                             % by CO2/brine interface area
                                             % in cell
    
      small = 3e-3;
      smooth_to_zero = @(x) value(x<small) .* ((x/small).^2 .* (- 2 * x / small + 3)) + value(x>=small) * 1;
      
      s_fac = smooth_to_zero(sG); % approximately one, but goes to 0 for very small values of sG
      rs_eps = (rsSat - rs) ./ rsSat;
      rs_fac = smooth_to_zero(rs_eps);
      rate = rate .* s_fac .* rs_fac; % smoothly turn down dissolution rate when
                                      % sG goes to 0, or when dissolved value
                                      % approaches maximum.
    
      % Conservation equation for dissolved CO2
      eqs{1} = (op.pv/dt) .* (pvMult  .* rs  .* bW  .* sW - pvMult0 .* rs0 .* bW0 .* sW0) + ...
               op.Div(rsbWvW) - rate;
    
      % accounting for dissolved gas in brine extracted from producer wells
      cqsw(cqsw > 0) = 0; % only take into account producer wells (inflow
                          % assumed not to contain any disgas)
      eqs{1}(wc) = eqs{1}(wc) - cqsw;
    
      % Taking boundary conditions and sources into account 
      
      flds = fieldnames(src);
      [vals, cells] = deal(cell(1, numel(flds)));
      for i = 1:numel(flds)
          obj = src.(flds{i});
          bcCells = obj.sourceCells;
          if ~isempty(bcCells)
              eqs{1}(bcCells) = eqs{1}(bcCells) - rs(bcCells).*obj.phaseMass{1}./rhoW(bcCells);
          end
          vals{i} = rs(bcCells).*obj.phaseMass{1}./rhoW(bcCells);
          cells{i} = bcCells;
      end
    
      % Modify conservation equation to ensure dissolved CO2 remains within
      % valid boundaries
    
      % Computing minimum allowed residual saturation per cell (the upper region where
      % both CO2 and brine phase is present should have saturated brine)
      min_rs = minRsHybrid(p, sG, sGmax, f, G, state);
    
      % Computing a hypothetic equation involving only minimum allowed
      % saturation.  We use this equation to identify where the real solution
      % would go below minimal allowed saturation.  Lock saturation of the
      % corresponding cells to the minimum value (but not lower).
      tmp = (op.pv / dt) .* ( value(pvMult) .* (value(min_rs) .* value(bW) ) - ...
                             pvMult0 .* (rs0 .* bW0 .* sW0 ) ) + ...
            op.Div(value(rsbWvW)) - value(rate);
      tmp(wc) = tmp(wc) - value(cqsw);
      for i = 1:numel(vals)
          tmp(cells{i}) = tmp(cells{i}) - value(vals{i});
      end
    
      ilow = tmp > -sqrt(eps); % If so, then min_rs is larger than the
                               % solution of eqs{1}, in other words, the
                               % solution of eqs{1} is smaller than the
                               % allowed value.  We have to modify eqs{1}.
      if any(ilow)
         % force value of 'rs' in these cells to equal 'min_rs'
         eqs{1}(ilow) = rs(ilow) .* sW(ilow) - min_rs(ilow);
         eqs{1}(ilow) = eqs{1}(ilow) .* op.pv(ilow)/dt;
      end
    
      % Identify cells with fully-saturated brine, and ensure saturation does
      % not rise further
      is_sat = (value(rs) >= value(rsSat)) & (eqs{1} < 0);
      eqs{1}(is_sat) = (rs(is_sat) - rsSat(is_sat)) .* op.pv(is_sat)/dt;
    
    
      % Equation for changes to sGmax
    
      % Default equation: force 'sGmax' to equal 'sG'
      eqs{2} = (sGmax - sG) .* (op.pv / dt);
    
      % In the following two special cases, sGmax will end up being higher
      % than sG, since any residual saturation is not fully eaten away by
      % dissolution (either due to rate being too slow, or due to saturation
      % of CO2 in brine being reached).
      swr = f.krPts.w(1);
      snr = f.krPts.g(1);
      tmp = (op.pv / dt) .* ...
            pvMult .* bG .* (sGmax - sGmax0) * ...
            snr ./ (1-swr);
      tmp2 = (op.pv / dt) .* ...
             pvMult .* value(bG) .* (value(sG) - sGmax0) * ...
             snr ./ (1 - swr) + rate;
    
      % Special case 1: New state remains UNSATURATED, but dissolution too slow to
      % deplete all residual saturation
      ix = (tmp2 < 0) & ~is_sat ;
      if any(ix)
         eqs{2}(ix) = tmp(ix) + rate(ix);
      end
    
      % Special case 2: New state reaches SATURATED value before all residual
      % saturation has been depleted
      ix = is_sat & (sGmax < sGmax0) & (sGmax > sG);
      if any(ix)
         eqs{2}(ix) = tmp(ix) + (rs(ix) - rs0(ix));
      end


      % Instantaneous dissolution model   
      isFine = ~isVE;
      %isFine = 1:model.G.cells.num;
%       eqs{1}(isFine) = sG(isFine) + 0*sG0(isFine);
%       eqs{1}(isSat & isFine) = rs(isSat & isFine) - rsSat(isSat & isFine); % set dissolved saturation rs to maximum possible (rsSat) for all cells with any CO2 (sG > 0) -> this CO2 will dissolve to the maximum possible instantaneously
%       eqs{2}(isFine)        = sGmax(isFine) - max(sGmax0(isFine), sG(isFine)); % standard equation for sMax -> no dependency on time-derivative (as is the case for rate-driven equations) since changes happen instantaneously
%       
end




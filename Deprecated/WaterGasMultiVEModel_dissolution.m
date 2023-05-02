classdef WaterGasMultiVEModel_dissolution < TwoPhaseWaterGasModelDissolution
    % Multi-VE version of the gas-water model from the co2lab module
    properties        
        sealingFaces
        sealingCells
        pe_rest
    end
    
    methods
        
        function model = WaterGasMultiVEModel_dissolution(G, rock, fluid, varargin)
            opt = struct('sealingFaces', [], 'sealingCells', [], 'pe_rest', 0);
            opt = merge_options(opt, varargin{:});
          
            % Modify getEquations to use
            % equationsWaterGasMultiVE_dissolution.
            % Then modify some of the functions to be applicable for VE cells
            model = model@TwoPhaseWaterGasModelDissolution(G, rock, fluid, varargin{:});                      
                      
            model.G.parent.sealingCells = opt.sealingCells; 
            model.G.sealingCells = model.G.partition(opt.sealingCells);
                        
            model.fluid.hys = ~isempty(regexp(func2str(model.fluid.krG), 'Hysteresis', 'match'));
            model.fluid.pe_rest = opt.pe_rest;
            
            model = model.setupOperators(); % recompute operators in case assignment to sealing cells did modify some operators 
            
            model.useCNVConvergence = true;                                  
            model = addCoarseOperatorsMultiVE(model);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                drivingForces, ...
                varargin)
            [problem, state] = equationsWaterGasMultiVE_dissolution(model, state0, state , dt , ...
                drivingForces, ...
                varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
           [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                        drivingForces);

           sg          = state.s(:,2);
           sg          = min(1-model.fluid.krPts.w(1), max(0, sg)); %(1-model.fluid.res_water, max(0,sg)); @@
           state.s     = [1-sg, sg];    
           
           state.sGmax = min(1-model.fluid.krPts.w(1), state.sGmax);
           state.sGmax = max(0,state.sGmax);
           state.sGmax = max(state.sGmax,sg); % Won't this overwrite the new governing (closure) equation for sGmax for dissolution model?

           % Specific update of rs for hybrid model:
           %isVE = model.G.cells.discretization > 1;
           isVE = true(model.G.cells.num, 1);
            if isfield(model.fluid, 'dis_rate')     
                f = model.fluid;
                % Rate-driven dissolution for VE cells
                 min_rs      = minRsHybrid(state.pressure, state.s(:,2), state.sGmax, ...
                                            f, model.G, state, isVE);
                 min_rs      = min_rs./state.s(isVE,1);
                 state.rs(isVE)    = max(min_rs,state.rs(isVE));
                 state.rs(isVE)    = min(state.rs(isVE), f.rsSat(state.pressure(isVE)));         
             
                % Instantaneous dissolution
                 diff = 1e-3; % @@  necessary for convergence in some cases
                 state.rs(~isVE) = min(state.rs(~isVE), model.fluid.rsSat(state.pressure(~isVE)) + diff);             
            end

        end
       
    end
end

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

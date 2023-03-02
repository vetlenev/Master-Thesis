classdef WaterGasMultiVEModel_PEBI < TwoPhaseWaterGasModel_PEBI
    % Multi-VE version of the gas-water model from the co2lab module
    properties
        
    end
    
    methods
        
        function model = WaterGasMultiVEModel_PEBI(G, rock, fluid, varargin)
            opt = struct('sealingFaces', [], 'sealingCells', [], 'pe_rest', 0);
            opt = merge_options(opt, varargin{:});
                       
            % Update to be dependent on TwoPhaseWaterGasModelHys instead,
            % since this is tailored for dissolution.
            % Modify getEquations to use
            % equationsWaterGasMultiVE_dissolution.
            % Then modify some of the functions to be applicable for VE cells
            model = model@TwoPhaseWaterGasModel_PEBI(G, rock, fluid, varargin{:});                      
                      
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
            [problem, state] = equationsWaterGasMultiVE_test(model, state0, state , dt , ...
                drivingForces, ...
                varargin{:});
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@TwoPhaseWaterGasModel_PEBI(model, state, problem, dx, ...
                                                drivingForces);

            if isfield(model.fluid, 'dis_rate')
              % The model includes dissolution
              if model.fluid.dis_rate > 0
                 % rate-driven dissolution
                 f           = model.fluid;
                 min_rs      = minRs(state.pressure,state.s(:,2),state.sGmax,f,model.G);
                 min_rs      = min_rs./state.s(:,1);
                 state.rs    = max(min_rs,state.rs);
                 state.rs    = min(state.rs,f.rsSat(state.pressure));         
              else
                 % instantaneous dissolution
                 diff = 1e-3; % @@  necessary for convergence in some cases
                 state.rs = min(state.rs, model.fluid.rsSat(state.pressure) + diff);
              end
            end
            
        end

        function g = getGravityVector(model)
        % Get the gravity vector used to instantiate the model
        %
        % SYNOPSIS:
        %   g = model.getGravityVector();
        %
        %
        % PARAMETERS:
        %   model - Class instance
        %
        % RETURNS:
        %   g     - `model.G.griddim` long vector representing the gravity
        %           accleration constant along increasing depth.
        %
        % SEE ALSO:
        %   `gravity`
        g = getGravityVector@ReservoirModel(model);
        if isfield(model.G, 'griddim')                       
            %dims = 1:model.G.griddim;
            dims = 1:3;
        else
            dims = ':';
        end
        g = model.gravity(dims);
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

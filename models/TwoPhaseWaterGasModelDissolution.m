classdef TwoPhaseWaterGasModelDissolution < ReservoirModel
    properties
        disgas
    end
    % Two-phase gas and water model        
    % ============================================================================
    methods
        % ------------------------------------------------------------------------
        function model = TwoPhaseWaterGasModelDissolution(G, rock, fluid, varargin)
           
            opt = struct;
            %model = model@TwoPhaseWaterGasModel(G, rock, fluid, tsurf, tgrad, varargin{:}); 
            [opt, unparsed] = merge_options(opt, varargin{:});%#ok
   
              model@ReservoirModel(G, varargin{:});
              model.rock    = rock;
              model.fluid   = fluid;
              model.water   = true;
              model.gas     = true;
              model.oil     = false;
              model.gravity = gravity;
              model.disgas = true;

              %model.equation = @equationsWGdisgas;

              model = model.setupOperators(G, rock);
        end
    end
    
    methods

        % ------------------------------------------------------------------------
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            [problem, state] = equationsWGdisgas(model, state0, state , dt , ...
                                                 drivingForces             , ...
                                                 varargin{:});
        end

        function [fn, index] = getVariableField(model, name, varargin)
      
          switch(lower(name))
            case {'sgmax'}
              index = 1;
              fn = 'sGmax';
            case {'rs'}
              index = 1;
              fn = 'rs';
            otherwise
              [fn, index] = getVariableField@ReservoirModel(model, name, varargin{:});
          end

        end

% Think this is the same as default setupOperators provided by ReservoirModel        
%         function model = setupOperators(model, G, rock, varargin)        
%           T             = computeTrans(G, rock); 
%           cf            = G.cells.faces(:, 1); 
%           nf            = G.faces.num; 
%           T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);       
%     
%           pv = poreVolume(G, rock);         
%     
%           model.operators = setupOperatorsTPFA(G, rock, 'porv', pv, 'trans', T);
%           model.operators.T_all = T;
%        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
           [state, report] = updateState@ReservoirModel(model, state, problem, dx, ...
                                                        drivingForces);

           sg          = state.s(:,2);
           sg          = min(1, max(0, sg)); %(1-model.fluid.res_water, max(0,sg)); @@
           state.s     = [1-sg, sg];    
           
           state.sGmax = min(1,state.sGmax);
           state.sGmax = max(0,state.sGmax);
           %state.sGmax = max(state.sGmax,sg);

            if isfield(model.fluid, 'dis_rate')
              % The model includes dissolution
              if model.fluid.dis_rate > 0
                 % rate-driven dissolution
                 f           = model.fluid;
                 min_rs      = minRsHybrid(state.pressure, state.s(:,2), state.sGmax, ...
                                            f, model.G, state);
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
        
    end
    
end
% ============================================================================

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


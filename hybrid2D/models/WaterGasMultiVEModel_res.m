classdef WaterGasMultiVEModel_res < TwoPhaseWaterGasModel
    % Multi-VE version of the gas-water model from the co2lab module
    properties
        
    end
    
    methods
        
        function model = WaterGasMultiVEModel_res(G, rock, fluid, varargin)
            opt = struct('sealingFaces', [], 'sealingCells', []);
            opt = merge_options(opt, varargin{:});
                       
            model = model@TwoPhaseWaterGasModel(G, rock, fluid, nan, nan, varargin{:});                      
            %model.G.parent.sealingFaces = opt.sealingFaces;           
            model.G.parent.sealingCells = opt.sealingCells;           
            %model.G.sealingFaces = model.G.faces.fconn(model.G.faces.connPos(opt.sealingFaces):model.G.faces.connPos(opt.sealingFaces+1)-1);
            model.G.sealingCells = model.G.partition(opt.sealingCells);
            
            model = model.setupOperators(); % recompute operators in case assignment to sealing cells did modify some operators 
            
            model.useCNVConvergence = true;                                  
            model = addCoarseOperatorsMultiVE(model);
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, ...
                drivingForces, ...
                varargin)
            [problem, state] = equationsWaterGasMultiVE_res(model, state0, state , dt , ...
                drivingForces, ...
                varargin{:});
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

function reports = makeFineReports(Gt_all, Gs_all, cmaps, fmaps, ...
                                    G, states, ...
                                    rock, fluid, schedule, ...
                                    residual, trap_s, dh)
% This function does intermediate processing of simulation data in order to
% generate inventory plots using 'plotTrappingDistribution'.
% 
% Currently, only rate controlled wells are supported (not pressure-controlled).
% 
% SYNOPSIS:
%   function reports = makeFineReports(Gt, states, rock, fluid, schedule, residual, traps, dh)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt_all    - all top-surface grids located at bottom of sealing layers
%   Gs_all    - all subgrids bounded by sealing layers (semi-perm or
%               caprock)
%   cmaps     - cell array of mappings from subgrid cells to global cells
%   fmaps     - cell array of mappings from subgrid faces to global faces
%   G        - FINE grid with all semi-perm layers in Gti
%   states   - result from a simulation (cell array of states, including
%              initial state)
%   rock     - rock object used in the simulation
%   fluid    - fluid object used in the simulation
%   schedule - schedule used in the simulation (NB: only rate controlled
%              wells supported at present)
%   residual - residual saturations, on the form [s_water, s_co2]
%   traps    - trapping structure (from trapAnalysis of Gt)
%   trapsi    - cell array of trapping structures (from trapAnalysis of Gti)
%   dh       - subscale trapping capacity (empty, or one value per grid cell
%              of Gt)
%
% RETURNS:
%   reports - a structure array of 'reports', that can be provided to the
%   'plotTrappingDistribution' function.
%
% SEE ALSO:
%   `plotTrappingDistribution`.

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
   
   assert( numel(states) == numel(schedule.step.val)+1 , ...
       'Ensure the initial state has been included in the varargin ''states''.')
   sw = residual(1);
   sn = residual(2);
 
   % Get global unique cells for VE subgrids under semi-perm layers   
    cs_all = {};
    for s=1:numel(Gs_all)
       cs_all = cat(1, cs_all, cmaps{s});
    end        
    cmaps_all = vertcat(cmaps{:});
  
    [~, ~, kk] = gridLogicalIndices(G);
  
    max_massVE = zeros(numel(Gs_all), 1); % max mass obtained in each VE top-surface subgrid    
    
   for i = 1:numel(states)      
      state = states{i};
      
      if i == 1          
         % handle initial state manually
         S = state.s(:,2);
         S_max = S; % at first step, max sat is just current sat
         %h     = G.cells.height .* ((S .* (1-sw) - (S_max .* sn)) ./  ...
         %                          ((1-sw) .* (1 - sw - sn)));
         %h_max = G.cells.height .* (S_max ./(1-sw));                  
         
         ntg = ones(G.cells.num,1);
         if isfield(rock,'ntg')
            ntg = rock.ntg;
         end
         % Set initial total "injected" volume to be the initial saturation
         % distribution in domain (necessary to get correct leakage)
         tot_inj = G.cells.volumes .* rock.poro .* ntg ...
                        .* S .* fluid.rhoGS; % multiply by S instead of h since saturation is uniquely defined in each fine cell
         tot_inj = sum(tot_inj); %
         reports(i).t         = 0; %#ok
         reports(i).W         = []; %#ok
      else          
         %h = state.h;
         %h_max = state.h_T; % only top residual plume, remaining residual plumes handled separately                                        
          
         reports(i).t = sum(schedule.step.val(1:i-1)); %#ok
         reports(i).W = schedule.control(schedule.step.control(i-1)).W; %#ok
         
         assert(all(cellfun(@(x) strcmpi(x, 'rate'), {reports(i).W.type})));
         tot_inj = tot_inj + (sum([reports(i).W.val]) * schedule.step.val(i-1) * fluid.rhoGS) ...
                                * reports(i).W.status; % only add injection rate if well is turned on
      end
      
      reports(i).sol       = state; % fine states
      %reports(i).sol.h     = h; % hybrid mobile heights
      %reports(i).sol.h_max = h_max; % hybrid residual heights
      
      rs = 0;
      if isfield(state, 'rs')
         rs = state.rs;
      end
     
      % FIRST, we compute CO2 mass for entire full-dim grid (no removal of
      % cells). No categorization, just sum up net mass.     
      reports(i).G.masses = state.s(:,2) .* G.cells.volumes .* rock.poro .* fluid.rhoGS;
                   
      % SECOND, compute masses for VE subgrids under top surface
      for s=1:numel(Gt_all)
          Gt = Gt_all{s}; % top surface of subgrid for semi-perm layer
          Gs = Gs_all{s}; % full subgrid for semi-perm layer
          c_sub = cmaps{s};
          f_sub = fmaps{s};
          trap_i = trap_s{s};
          Gts = sprintf('Gts%d', s); 
          
          if i==numel(states)
             test = 0; 
          end
                
          reports(i).(Gts).masses    = massTrappingDistributionFineADI(Gt, Gs, G, ...
                                                        c_sub,  ...
                                                        reports(i).sol.pressure , ...
                                                        reports(i).sol.s(:,2)   , ...
                                                        reports(i).sol.s(:,1)   , ...                                                        
                                                        rock                    , ...
                                                        fluid                   , ...
                                                        trap_i                  , ...
                                                        dh                      , ...
                                                        'rs', rs); %#ok
                                                    
          mass_i = reports(i).(Gts).masses; % individual masses for this subgrid
          max_massVE(s) = max(max_massVE(s), sum(mass_i)); % keep track of max mass reached in this subgrid
          if s == 1
              mass_sub = zeros(size(mass_i)); % separate masses to be summed over subgrids
          end
          mass_sub = mass_sub + mass_i; % individual masses summed over subgrids
          
          leaked_i = max_massVE(s) - sum(mass_i); % leaked out of each subgrid separately
          reports(i).(Gts).masses = [reports(i).(Gts).masses, leaked_i]; % 0 -> leaked_i
      end
           
      cfh_rem = setdiff(G.cells.indexMap, cmaps_all);

      % Compute mass for remaining regions (inside cell-represented
      % semi-perm layers)
      mass_rem = massTrappingFine_Other(G, cfh_rem, reports(i).sol.pressure, ...
                                                        reports(i).sol.s(:,2), ...
                                                        reports(i).sol.s(:,1), ...
                                                        rock, fluid, 'rs', rs);   
      
      
      % Sum up masses BEFORE adding leaked mass
      reports(i).masses =  mass_sub + mass_rem;
                
      %leaked = tot_inj - sum(reports(i).Gt.masses) - tot_mass_subgrids;
      leaked = tot_inj - sum(reports(i).G.masses);
      reports(i).masses = [reports(i).masses, leaked]; % 0 -> leaked
      % ---------------------------------------------------
   end
   
end
% ----------------------------------------------------------------------------

function [h, h_max] = compute_plume_height(Gt, state, sw, sr)
    
    if isfield(state, 'sGmax')
       smax = state.sGmax; % we operate with dissolution
    else
       smax = state.smax(:,2); % no dissolution.  Current max = historical max
    end
    [h, h_max] = upscaledSat2height(state.s(:,2), smax, Gt, 'resSat', [sw, sr]);
end

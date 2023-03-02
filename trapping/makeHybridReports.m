function reports = makeHybridReports(Gts_all, Gs_all, cmaps, fmaps, ...
                                    fine_rem, ve_rem, ...
                                    Gh, states, ...
                                    rock, fluid, schedule, ...
                                    residual, trap_s, dh)
% This function does intermediate processing of simulation data in order to
% generate inventory plots using 'plotTrappingDistribution'.
% 
% Currently, only rate controlled wells are supported (not pressure-controlled).
% 
% SYNOPSIS:
%   function reports = makeReports(Gt, states, rock, fluid, schedule, residual, traps, dh)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt       - top surface grid of full formation (i.e. confining caprock)
%   Gti      - cell array of top surface grids of interior lowperm layers
%   Gh       - hybrid grid with all semi-perm layers in Gti
%   cmaps     - cell array of mappings from subgrid cells to global cells
%   fmap     - cell array of mappings from subgrid faces to global faces
%   nmap     - cell array of mappings from subgrid nodes to global nodes
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
  
    [~, ~, kk] = gridLogicalIndices(Gh.parent);   
    ph = Gh.partition;        
  
    max_massVE = zeros(numel(Gts_all), 1); % max mass obtained in each VE top-surface subgrid    
    
   for i = 1:numel(states)

      %[h, h_max] = compute_plume_height(Gt, states{i}, residual(1), residual(2));      
      state = states{i};
      
      if i == 1          
         % handle initial state manually
         S = states{i}.s(:,2);
         S_max = S; % at first step, max sat is just current sat
         h     = Gh.cells.height .* ((S .* (1-sw) - (S_max .* sn)) ./  ...
                                   ((1-sw) .* (1 - sw - sn)));
         h_max = Gh.cells.height .* (S_max ./(1-sw));
         
         h_B = inf;
         h_H = inf; Hh = inf;
         h_BH = inf; Hbh = inf;
         cBottom = [];
         cHorz = [];
         cBottomHorz = [];
         
         ntg = ones(Gh.cells.num,1);
         if isfield(rock,'ntg')
            ntg = rock.ntg;
         end
         % Set initial total "injected" volume to be the initial saturation
         % distribution in domain (necessary to get correct leakage)
         %tot_inj = Gh.cells.volumes .* rock.poro .* ntg .* (1-sw) .* h .* fluid.rhoGS; % (1-sw) to account for residual water
         %tot_inj = Gt.cells.volumes .* rock.poro .* ntg .* (1-sw) .* h .* fluid.rhoGS;         
         % --- Think this is OK: ---
         tot_inj = Gh.cells.volumes .* rock.poro .* ntg .* S .* fluid.rhoGS;              
         % -------------------------
         tot_inj = sum(tot_inj); %
         reports(i).t         = 0; %#ok
         reports(i).W         = []; %#ok
      else          
          h = state.h;
          h_max = state.h_T; % only top residual plume, remaining residual plumes handled separately                    
          
          cBottom = state.cBottom;%(ismember(state.cBottom, cth)); % important to choos only bottom VE cells part of subgrid Gt
          cHorz = state.cHorz;%(ismember(state.cHorz, cth));
          cBottomHorz = state.cBottomHorz;%(ismember(state.cBottomHorz, cth));
          
          h_B = state.h_B(cBottom); % only choose cells where flux actually migrates from bottom (h_B contains ALL bottom 
          h_H = state.hHi; % for cells where CO2 migrates from veHorizontal layers
          h_BH = state.hBHi; % for cells where CO2 migrates from combined veHorizontal / veBottom layers                   
          
          Hh = state.Hi;
          Hbh = state.BHi;
          
         reports(i).t = sum(schedule.step.val(1:i-1)); %#ok
         reports(i).W = schedule.control(schedule.step.control(i-1)).W; %#ok
         
         assert(all(cellfun(@(x) strcmpi(x, 'rate'), {reports(i).W.type})));
         tot_inj = tot_inj + (sum([reports(i).W.val]) * schedule.step.val(i-1) * fluid.rhoGS) ...
                                * reports(i).W.status; % only add injection rate if well is turned on
      end
      
      reports(i).sol       = states{i}; % hybrid states
      reports(i).sol.h     = h; % hybrid mobile heights
      reports(i).sol.h_max = h_max; % hybrid residual heights
      reports(i).sol.h_B   = h_B;
      reports(i).sol.h_H   = h_H;
      reports(i).sol.h_BH  = h_BH;
      reports(i).sol.Hh    = Hh;
      reports(i).sol.Hbh   = Hbh;
      
      rs = 0;
      if isfield(states{i}, 'rs')
         rs = states{i}.rs;
      end
      
      % Masses for hybrid gid computed in 5 steps:
      
      % FIRST, we compute CO2 mass for entire hybrid grid (no removal of
      % cells). No categorization, just sum up net mass.     
      reports(i).Gh.masses = states{i}.s(:,2) .* Gh.cells.volumes .* rock.poro .* fluid.rhoGS;
                    
      if i == numel(states)
         test = 0; 
      end
      % SECOND, compute masses for VE subgrids under top surface
      for s=1:numel(Gts_all)
          Gs = Gts_all{s}; % top surface grid for semi-perm layer s
          c_sub = cmaps{s};
          f_sub = fmaps{s};
          trap_i = trap_s{s};
          Gts = sprintf('Gts%d', s);                      
                
          reports(i).(Gts).masses    = massTrappingDistributionHybridADI(Gs, Gh, ...
                                                        c_sub,  ...
                                                        reports(i).sol.pressure , ...
                                                        reports(i).sol.s(:,2)   , ...
                                                        reports(i).sol.s(:,1)   , ...
                                                        reports(i).sol.h        , ...
                                                        reports(i).sol.h_max    , ...
                                                        reports(i).sol.h_B     , ...
                                                        reports(i).sol.h_H     , ...
                                                        reports(i).sol.h_BH    , ...
                                                        reports(i).sol.Hh      , ...
                                                        reports(i).sol.Hbh     , ...
                                                        cBottom                  , ...
                                                        cHorz                    , ...
                                                        cBottomHorz              , ...
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
           
      %cfh_rem = ph(fine_rem);
      fine_rem_hybrid = ph(fine_rem);
      ve_rem_hybrid = ph(ve_rem);  

      if ~isempty(fine_rem_hybrid)
            mass_fine = massTrappingFine_Other(Gh, fine_rem_hybrid, reports(i).sol.pressure, ...
                                                        reports(i).sol.s(:,2), ...
                                                        reports(i).sol.s(:,1), ...
                                                        rock, fluid, 'rs', rs);   
      else
          mass_fine = 0;
      end

      if i == 500
          test = 0;
      end
    
      if ~isempty(ve_rem_hybrid)
            mass_ve = massTrappingVE_Other(Gh, ve_rem, ve_rem_hybrid, reports(i).sol.pressure, ...
                                                        reports(i).sol.s(:,2), ...
                                                        reports(i).sol.s(:,1), ...
                                                        reports(i).sol.h, ...
                                                        reports(i).sol.h_max, ...
                                                        reports(i).sol.h_B     , ...
                                                        reports(i).sol.h_H     , ...
                                                        reports(i).sol.h_BH    , ...
                                                        reports(i).sol.Hh      , ...
                                                        reports(i).sol.Hbh     , ...
                                                        cBottom                  , ...
                                                        cHorz                    , ...
                                                        cBottomHorz              , ...
                                                        rock, fluid, 'rs', rs); 
%             mass_ve = massTrappingVE_OtherNew(Gh, ve_rem_hybrid, reports(i).sol.pressure, ...
%                                                             reports(i).sol.s(:,2), ...
%                                                             reports(i).sol.s(:,1), ...
%                                                             rock, fluid, 'rs', rs);
      else
          mass_ve = 0;
      end     
      
      % Sum up masses BEFORE adding leaked mass
      reports(i).masses =  mass_sub + mass_fine + mass_ve;
          
      leaked = tot_inj - sum(reports(i).Gh.masses);
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

function reports = makeHybridReports(Gts_all, cmaps, fmaps, ...
                                    Gtf_all, cmapf, fmapf, ...
                                    fine_rem, ...
                                    Gh, states, ...
                                    rock, fluid, schedule, ...
                                    residual, trap_s, trap_f, dh)
% NB: Change cmaps to cmaps
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
 
   % Get global unique cells for VE subgrids under semi-perm   
    cs_all = {};
    for s=1:numel(Gts_all)
       cs_all = cat(1, cs_all, cmaps{s});
    end
    
    %cmap_glob = setdiff(Gt.parent.cells.indexMap, vertcat(cs_all{:})); % extract all cells from 3D grid associated uniquely with Gt (i.e. cells not part of any Gti)    
    [~, ~, kk] = gridLogicalIndices(Gh.parent);
    
    ph = Gh.partition;
    % Get topsurface-hybrid mapping
    % --- Alternative 1 ---
%     kkt = find(kk == min(kk));
%     ctf_glob = cmap_glob(ismember(cmap_glob, kkt)); % VE regions of global top surface
%     cth_glob = ph(ctf_glob); % same but hybrid indices 
   
    max_massVE = zeros(numel(Gts_all), 1); % max mass obtained in each VE top-surface subgrid
    max_massFine = zeros(numel(Gtf_all), 1); % max mass obtained in each fine top-surface subgrid
    
   for i = 1:numel(states)

      %[h, h_max] = compute_plume_height(Gt, states{i}, residual(1), residual(2));      
      state = states{i};
      
      if i == 1          
         % initial state can contain non-zero co2 saturations (or heights)
         S = states{i}.s(:,2);
         %S = S(cth);
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
         tot_inj = Gh.cells.volumes .* rock.poro .* ntg .* (1-sw) ...
                        .* h .* fluid.rhoGS;
         tot_inj = sum(tot_inj); %
         reports(i).t         = 0; %#ok
         reports(i).W         = []; %#ok
      else          
          h = states{i}.h;
          h_max = state.h_T; % only top residual plume, remaining residual plumes handled separately                    
          
          cBottom = state.cBottom;%(ismember(state.cBottom, cth)); % important to choos only bottom VE cells part of subgrid Gt
          cHorz = state.cHorz;%(ismember(state.cHorz, cth));
          cBottomHorz = state.cBottomHorz;%(ismember(state.cBottomHorz, cth));
          
          h_B = state.h_B(cBottom); % only choose cells where flux actually migrates from bottom
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
      
      % SECOND, we only compute CO2 mass for global top surface, using same top surface
      % but removing cells that are part of other top surface subgrid and other fine cells.
%       reports(i).Gt.masses    = massTrappingDistributionHybridADI(Gt, Gh, ...
%                                                         cth_glob,  ...
%                                                         reports(i).Gt.sol.pressure , ...
%                                                         reports(i).Gt.sol.s(:,2)   , ...
%                                                         reports(i).Gt.sol.s(:,1)   , ...
%                                                         reports(i).Gt.sol.h        , ...
%                                                         reports(i).Gt.sol.h_max    , ...
%                                                         reports(i).Gt.sol.h_B     , ...
%                                                         reports(i).Gt.sol.h_H     , ...
%                                                         reports(i).Gt.sol.h_BH    , ...
%                                                         reports(i).Gt.sol.Hh      , ...
%                                                         reports(i).Gt.sol.Hbh     , ...
%                                                         cBottom                  , ...
%                                                         cHorz                    , ...
%                                                         cBottomHorz              , ...
%                                                         rock                    , ...
%                                                         fluid                   , ...
%                                                         trap_s                   , ...
%                                                         dh                      , ...
%                                                         'rs', rs); %#ok                     
%       
%       mass_top_glob = reports(i).Gt.masses; % individual masses for each subgrid
%       max_massGlob = max(max_massGlob, sum(mass_top_glob)); % keep track of max mass reached in this subgrid     
%           
%       leaked_top = max_massGlob - sum(mass_top_glob); % leaked out of each subgrid separately
%       reports(i).Gt.masses = [reports(i).Gt.masses, leaked_top];
                                                    
      % --- Loop through interior semi-permeable layers ---         
      % THIRD, compute masses for VE subgrids under top surface
      for s=1:numel(Gts_all)
          Gs = Gts_all{s}; % top surface subgrid for semi-perm layer
          cs = cmaps{s};
          fs = fmaps{s};
          trap_i = trap_s{s};
          Gts = sprintf('Gts%d', s);
                                     
          kkf = kk(cs);       
          kkf = cs(kkf == min(kkf));          
          csf = cs(ismember(cs, kkf));         
          csh = ph(csf); % for new extraction of subgrid there will be no overlap for top surface subgids        
                
          reports(i).(Gts).masses    = massTrappingDistributionHybridADI(Gs, Gh, ...
                                                        csh,  ...
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
                                                        trap_i                   , ...
                                                        dh                      , ...
                                                        'rs', rs); %#ok
                                                    
          mass_i = reports(i).(Gts).masses; % individual masses for this subgrid
          max_massVE(s) = max(max_massVE(s), sum(mass_i)); % keep track of max mass reached in this subgrid
          if s == 1
              mass_subVE = zeros(size(mass_i)); % separate masses to be summed over subgrids
          end
          mass_subVE = mass_subVE + mass_i; % individual masses summed over subgrids
          
          leaked_i = max_massVE(s) - sum(mass_i); % leaked out of each subgrid separately
          reports(i).(Gts).masses = [reports(i).(Gts).masses, leaked_i]; % 0 -> leaked_i
      end
           
      mass_subFine = zeros(size(mass_subVE));      
      % FOURTH, compute masses for fine regions under top surface subgrids
      for s=1:numel(Gtf_all)
          Gf = Gtf_all{s}; % top surface subgrid for semi-perm layer
          cf = cmapf{s};
          ff = fmapf{s};
          trap_i = trap_f{s};
          Gtf = sprintf('Gtf%d', s-1); % s-1 so global top surface gets 0-index
                                     
          kkf = kk(cf);       
          kkf = cf(kkf == min(kkf));          
          cff = cf(ismember(cf, kkf));         
          cfh = ph(cff); % for new extraction of subgrid there will be no overlap for top surface subgids        
                
          reports(i).(Gtf).masses    = massTrappingFine_TopSurface(Gh, cfh,  ...
                                                        reports(i).sol.pressure , ...
                                                        reports(i).sol.s(:,2)   , ...
                                                        reports(i).sol.s(:,1)   , ...                                                       
                                                        rock                    , ...
                                                        fluid                   , ...
                                                        trap_i                   , ...              , ...
                                                        'rs', rs); %#ok
                                                    
          mass_i = reports(i).(Gtf).masses;
          max_massFine(s) = max(max_massFine(s), sum(mass_i));
          mass_subFine = mass_subFine + mass_i;
          
          leaked_i = max_massFine(s) - sum(mass_i); % leaked out of each subgrid separately
          reports(i).(Gtf).masses = [reports(i).(Gtf).masses, leaked_i];
      end
      % ------------------
                                                    
      % FIFTH, compute masses in fine parts of grid, not accounted for in
      % calculations for top surface
      kkf = kk(fine_rem);       
      kkf = cf(kkf == min(kkf));          
      cff = cf(ismember(cf, kkf));         
      cfh_rem = ph(cff); % for new extraction of subgrid there will be no overlap for top surface subgids        

      mass_fine = massTrappingFine_Other(Gh, cfh_rem, reports(i).sol.pressure, ...
                                                        reports(i).sol.s(:,2), ...
                                                        reports(i).sol.s(:,1), ...
                                                        rock, fluid, 'rs', rs);   
      
      
      % Sum up masses BEFORE adding leaked mass
      reports(i).masses =  mass_subVE + ...
                            mass_subFine + ...
                            mass_fine;
                
      %leaked = tot_inj - sum(reports(i).Gt.masses) - tot_mass_subgrids;
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

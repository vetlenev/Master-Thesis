function reports = makeHybridReports(Gt, cmap_glob, fmap_glob, ...
                                    Gti, cmaps, fmaps, ...
                                    Gtf, cmapf, fmapf, ...
                                    fine_rem, ...
                                    Gh, states, ...
                                    rock, fluid, schedule, ...
                                    residual, traps, trapsi, dh)
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
    for s=1:numel(Gti)
       cs_all = cat(1, cs_all, cmaps{s});
    end
    
    ct_all = setdiff(Gt.parent.cells.indexMap, vertcat(cs_all{:})); % extract all cells from 3D grid associated uniquely with Gt (i.e. cells not part of any Gti)
    [~, ~, kk] = gridLogicalIndices(Gt.parent);
    
    ph = Gh.partition;
    % Get topsurface-hybrid mapping
    % --- Alternative 1 ---
    kkt = find(kk == min(kk));
    ctf = ct_all(ismember(ct_all, kkt)); % fine cells at global top surface
    cth = ph(ctf); % fine cells at global top surface with hybrid indices 
   
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
         
         ntg = ones(Gt.cells.num,1);
         if isfield(rock,'ntg')
            ntg = rock.ntg;
         end
         % Set initial total "injected" volume to be the initial saturation
         % distribution in domain (necessary to get correct leakage)
         tot_inj = Gt.cells.volumes .* rock.poro(cth) .* ntg .* (1-sw) ...
             .* h(cth) .* fluid.rhoGS;
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
      
      reports(i).Gt.sol       = states{i}; % hybrid states
      reports(i).Gt.sol.h     = h; % hybrid mobile heights
      reports(i).Gt.sol.h_max = h_max; % hybrid residual heights
      reports(i).Gt.sol.h_B   = h_B;
      reports(i).Gt.sol.h_H   = h_H;
      reports(i).Gt.sol.h_BH  = h_BH;
      reports(i).Gt.sol.Hh    = Hh;
      reports(i).Gt.sol.Hbh   = Hbh;
      
      rs = 0;
      if isfield(states{i}, 'rs')
         rs = states{i}.rs;
      end
      
      % First, we compute CO2 mass for entire hybrid grid (no removal of cells)
      reports(i).Gh.masses = states{i}.s(:,2) .* Gh.cells.volumes .* rock.poro .* fluid.rhoGS;
      
      % Next, we only compute CO2 mass for unique parts of top surface, using same top surface
      % but removining cells that are part of another top surface grid and other fine cells.
      reports(i).Gt.masses    = massTrappingDistributionHybridADI(Gt, Gh, ...
                                                        cth,  ...
                                                        reports(i).Gt.sol.pressure , ...
                                                        reports(i).Gt.sol.s(:,2)   , ...
                                                        reports(i).Gt.sol.s(:,1)   , ...
                                                        reports(i).Gt.sol.h        , ...
                                                        reports(i).Gt.sol.h_max    , ...
                                                        reports(i).Gt.sol.h_B     , ...
                                                        reports(i).Gt.sol.h_H     , ...
                                                        reports(i).Gt.sol.h_BH    , ...
                                                        reports(i).Gt.sol.Hh      , ...
                                                        reports(i).Gt.sol.Hbh     , ...
                                                        cBottom                  , ...
                                                        cHorz                    , ...
                                                        cBottomHorz              , ...
                                                        rock                    , ...
                                                        fluid                   , ...
                                                        traps                   , ...
                                                        dh                      , ...
                                                        'rs', rs); %#ok
                                                    
      % Then compute masses in fine parts of grid, not accounted for in
      % calculations for top surface
      fine_masses = fineMassTrappingHybrid(
      
      
      % --- Loop through interior semi-permeable layers ---
      tot_mass_subgrids = 0;
      mass_subgrids = zeros(size(reports(i).Gt.masses));
      max_mass = zeros(numel(Gti), 1); % max mass obtained in each subgrid
      
      for s=1:numel(Gti)
          Gs = Gti{s}; % top surface subgrid for semi-perm layer
          cs = cmaps{s};
          fs = fmaps{s};
          trap_i = trapsi{s};
          Gts = sprintf('Gt%d', s);
                    
          cs_other = cs_all;
          cs_other{s} = [];
          % Extract unique subcells and map from fine to hybrid indices          
          %cs_unique = setdiff(cs, vertcat(cs_other{:}));          
          kkf = kk(cs);       
          kkf = cs(kkf == min(kkf));          
          csf = cs(ismember(cs, kkf));         
          csh = ph(csf); % for new extraction of subgrid there will be no overlap for top surface subgids        
                
          reports(i).(Gts).masses    = massTrappingDistributionHybridADI(Gs, Gh, ...
                                                        csh,  ...
                                                        reports(i).Gt.sol.pressure , ...
                                                        reports(i).Gt.sol.s(:,2)   , ...
                                                        reports(i).Gt.sol.s(:,1)   , ...
                                                        reports(i).Gt.sol.h        , ...
                                                        reports(i).Gt.sol.h_max    , ...
                                                        reports(i).Gt.sol.h_B     , ...
                                                        reports(i).Gt.sol.h_H     , ...
                                                        reports(i).Gt.sol.h_BH    , ...
                                                        reports(i).Gt.sol.Hh      , ...
                                                        reports(i).Gt.sol.Hbh     , ...
                                                        cBottom                  , ...
                                                        cHorz                    , ...
                                                        cBottomHorz              , ...
                                                        rock                    , ...
                                                        fluid                   , ...
                                                        trap_i                   , ...
                                                        dh                      , ...
                                                        'rs', rs); %#ok
          %leaked = tot_inj - sum(reports(i).Gt.masses);
          %reports(i).Gt.masses = [reports(i).Gt.masses, leaked];
          % Only calculate leakage for interior layers at boundary of
          % parent grid
          boundary_faces_G = boundaryFaces(Gt.parent);
          at_boundary = ismember(fs, boundary_faces_G);
          if any(at_boundary)
             fs_boundary = fs(at_boundary);
             %leaked = leaked - 
          end
          
          % FIX THIS:
          %leaked = 0;
          %reports(i).(Gts).masses = [reports(i).(Gts).masses, leaked];
          
          mass_i = reports(i).(Gts).masses; % individual masses for each subgrid
          max_mass(s) = max(max_mass(s), sum(mass_i)); % keep track of max mass reached in this subgrid
          mass_subgrids = mass_subgrids + mass_i; % individual masses summed over subgrids
          tot_mass_subgrids = tot_mass_subgrids + sum(mass_i); % net CO2 mass summed over subgrids
          
          leaked_i = max_mass(s) - sum(mass_i); % leaked out of each subgrid separately
          reports(i).(Gts).masses = [reports(i).(Gts).masses, leaked_i]; % 0 -> leaked_i
      end
      
      % Sum up masses BEFORE adding leaked mass
      reports(i).masses = reports(i).Gt.masses + mass_subgrids;
                
      %leaked = tot_inj - sum(reports(i).Gt.masses) - tot_mass_subgrids;
      leaked = tot_inj - sum(reports(i).Gh.masses);      
      reports(i).Gt.masses = [reports(i).Gt.masses, leaked]; % 0 -> leaked
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

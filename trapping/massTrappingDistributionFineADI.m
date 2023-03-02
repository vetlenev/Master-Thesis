function [masses, masses0] = massTrappingDistributionFineADI(Gt, Gs, G, c_fine, p, sG, sW, ...
                                                                 rock, fluidADI, traps, dh, varargin)
% Compute the trapping status distribution of CO2 in full-dimensional
% model.
%
% SYNOPSIS:
%   masses = massTrappingDistributionVEADI(Gt, p, sW, sG, h, h_max, ...
%                              rock, fluidADI, sr, sw, trapstruct)
%   masses = massTrappingDistributionVEADI(Gt, p, sW, sG, h, h_max, ...
%                              rock, fluidADI, sr, sw, trapstruct, 'rs',rs)
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gt         - Extracted top-surface grid
%   ?Gti?        - cell array of semi-permeable top surface grids
%   G          - fine grid
%   cmap       - cell array of local-to-global cell mapping for 3D grid corresponding to each Gti
%   p          - pressure, one value per cell of grid
%   sW         - water saturation, one value per cell of HYBRID grid
%   sG         - gas saturation, one value per cell of grid
%   h          - gas height below top surface, one value per cell of grid
%   h_max      - maximum historical gas height, one value per cell of grid
%   rock       - rock parameters corresponding to 'Gt'
%   fluidADI   - ADI fluid object (used to get densities and compressibilities)
%   sr         - gas residual saturation (scalar)
%   sw         - liquid residual saturation (scalar)
%   trapstruct - trapping structure
%   dh         - subtrapping capacity (empty, or one value per grid cell of Gt)
%   varargin   - optional parameters/value pairs.  This currently only
%                includes the option 'rs', which specifies the amount of
%                dissolved CO2 (in its absence, dissolution is ignored).
%
% RETURNS:
%   masses - vector with 7 components, representing:
%            masses[1] : mass of dissolved gas, per cell
%            masses[2] : mass of gas that is both structurally and residually trapped
%            masses[3] : mass of gas that is residually (but not structurally) trapped
%            masses[4] : mass of non-trapped gas that will be residually trapped
%            masses[5] : mass of structurally trapped gas, not counting the gas that 
%                        will eventually be residually trapped
%            masses[6] : mass of subscale trapped gas (if 'dh' is nonempty)
%            masses[7] : mass of 'free' gas (i.e. not trapped in any way)
%   masses0 (optional) - masses in terms of one value per grid cell of Gt

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

    opt.rs = 0;
    opt    = merge_options(opt, varargin{:});        
       
    [ii, jj, kk] = gridLogicalIndices(G);
    % --- Collect all cells with similar ii and jj indices ---
    %accumarray([ii;jj], traps.z_spill_loc, [], @mode);
    fine_cells_sub = cell(Gt.cells.num,1);    
    fine_z = cell(Gt.cells.num,1);
    %z_spill_loc = cell(Gt.cells.num,1);
    if Gt.cells.num == 1
        Gt.columns.cells = Gt.columns.cells'; % transpose if only one top surface cell
    end
    for i=1:Gt.cells.num
        fine_cells_sub{i} = Gt.columns.cells(Gt.cells.columnPos(i):Gt.cells.columnPos(i+1)-1, :); % fine cells part of column associated with cell i from top-surface grid             
        fine_z{i} = Gs.cells.centroids(fine_cells_sub{i}, 3); % z-coord of fine cells in column from SUBGRID
        % OR USE TOP OF CELL AND NOT CENTROID, AS IN MAKE_HYBRID ?!?
    end
    fine_z_all = vertcat(fine_z{:});
    % --------------------------------------------------------
   
    %z = G.cells.centroids(:,3);
    %z = z(c_fine);           
    
    % Extracting relevant information from 'sol'
    sw=fluidADI.krPts.w(1);%liquid residual saturation (scalar)
    sr=fluidADI.krPts.g(1);%gas residual saturation (scalar)
    %SW = sW(c_fine);   
    %SG = sG(c_fine);
    
    rs = opt.rs;
    %p = p(c_fine);
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro .* G.cells.volumes .* pvMult; % since cells ae full-dim, use Gh.cells.volumes, not volumes from top surface
    if isfield(rock,'ntg')
       pv = rock.poro .* G.cells.volumes .* rock.ntg .* pvMult;
       %effective area accounting for possible net-to-gross data
    end
    rhoCO2   = fluidADI.rhoGS .* fluidADI.bG(p);
    gasPhase = sum(pv .* (rhoCO2 .* sG));        
    
    
    zt = cell(Gt.cells.num, 1);
    if isfield(traps, 'z_spill_loc')
        % trap heights -> only calculated from geometry, not h or h_max       
        for i=1:Gt.cells.num
            zt{i} = traps.z_spill_loc(i);
        end
    else
        for i=1:Gt.cells.num
            tr   = traps.trap_regions(i); % trap region for given cell of top-surface grid
            %tr_u = setdiff(unique(tr), 0);
            %tr_z = [min(fine_z_all); traps.trap_z(tr_u)]; % trap region zero -> spills out of domain -> can't become structurally trapped -> set spill region to minimum possible (min(all_fine_z)) to avoid assigning the cells as structurally trapped
            tr_z = [min(fine_z_all); traps.trap_z];
            tr   = 1 + tr; % shift so trap region 0 (no trap) extracts min(fine_z_all) as spill-point (which gives no leakage)
            %tr(tr>numel(tr_z)) = 1; % to ensure all trap regions with no calculated heights are assigned zero spill region height
            zt{i} = tr_z(tr); % should be scalar (same for all cells in column)
        end
    end     
    
    freeStruc = 0; resStruc = 0; freeRes = 0;
    freeMov = 0; resTrap = 0; resDis = 0;
    
    for i=1:Gt.cells.num % number of cells in top surface grid
        fine_z_global = c_fine(fine_cells_sub{i}); % map from subgrid indices to global indices
        %cols_i = cols_u(i);
        %col_idx = cols == cols_i;
        zi = fine_z{i};
        SG_i = sG(fine_z_global); % only select saturations from (global) fine cells in current trap region / column
        SW_i = sW(fine_z_global);
        p_i = p(fine_z_global);
        pv_i = pv(fine_z_global);
        rho_pv = fluidADI.rhoGS.*pv_i;
        
        strucTrapped_i = zi < zt{i};
        % -----
        res_buff = 1.5; % buffer factor for CO2 to be residually trapped (necessary since CO2 sat will converge towards sr at infinite time due to relperm going to zero)
        % -----
        freePlume_i = zi >= zt{i} & SG_i > res_buff*sr; % cells in mobile plume
        resPlume_i = zi >= zt{i} & SG_i <= res_buff*sr; % cells in immobilized plume   
        free = zi >= zt{i};
        
        freeStruc = freeStruc + sum((max(SG_i, sr) - sr).*rho_pv.* strucTrapped_i); % structural plume
        resStruc = resStruc + sum(min(SG_i, sr).*rho_pv .* strucTrapped_i); % structural residual
        freeRes = freeRes + sum(sr.*rho_pv .* freePlume_i); % residual in plume (we know SG_i > sr, so residual value can be fixed at sr)
        %freeRes = freeRes + sum(min(SG_i, sr).*rho_pv .* free);
        freeMov = freeMov + sum((SG_i - sr).*rho_pv .* freePlume_i); % free plume
        %freeMov = freeMov + sum(max(SG_i-sr,0).*rho_pv .* free);
        resTrap = resTrap + sum(SG_i.*rho_pv .* resPlume_i); % we know that SG_i <= sr, so no need to compute min(SG_i, sr)
        resDis = resDis + sum(rho_pv.* (rs .* fluidADI.bW(p_i) .* SW_i));        
    end
            
    %subtrap   = sum((hm_sub * sr + h_sub * (1 - sr - sw)) .* pv .* rhoCO2);
    subtrap = 0;
          
    masses    = max([value(resDis), value(resStruc), value(resTrap), ...
                     value(freeRes), value(freeStruc), value(subtrap), ...
                     value(freeMov)], 0); % may be ADI variables
         
    if(abs(sum(masses(2:end))-gasPhase) > 1e-3 * gasPhase)
        %disp('There is a mismatch between mass calculations');
        %fprintf('Mass: %d. Gas: %d\n', abs(sum(masses(2:end))-gasPhase), 1e-3 * gasPhase);       
    end
    
    if nargout > 1 % NB: NOT IMPLEMENTED YET !
        % values one per cell of grid Gt
        strucVol_0  = 0;
        plumeVol_0  = 0;
        resStruc_0  = 0;
        freeStruc_0 = 0;
        freeRes_0   = 0;
        freeMov_0   = 0;
        resTrap_0   = 0;
        resDis_0    = 0;
        subtrap_0   = 0;

        masses0 = [{resDis_0}, {resStruc_0}, {resTrap_0}, {freeRes_0}, ...
                    {freeStruc_0}, {subtrap_0}, {freeMov_0}]; % may be ADI vars
    end
end
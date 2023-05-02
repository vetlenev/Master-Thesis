function [masses, masses0] = massTrappingDistributionHybridADI(Gt, Gh, c_sub, p, sG, sW, h, h_max, ...
                                                                 h_B, h_H, h_BH, Hh, Hbh, cB, cH, cBH, ...
                                                                 rock, fluidADI, traps, dh, varargin)
% Compute the trapping status distribution of CO2 in each cell of global
% top-surface grid, removing cells part of top surface of semi-perm layers
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
%   Gt         - Top surface grid
%   ?Gti?      - cell array of semi-permeable top surface grids
%   Gh         - hybrid grid
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
%   masses_0 (optional) - masses in terms of one value per grid cell of Gt

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
    
    ph = Gh.partition;
    %c_subh = ph(c_sub);
    [~, ~, kk] = gridLogicalIndices(Gh.parent);
    
      ve_mask = Gh.cells.discretization(ph(c_sub)) > 1;          
      c_ve = c_sub(ve_mask); % ve cells of top surface subgrid
      c_fine = c_sub(~ve_mask); % fine cells of top surface subgrid

      fine_col = Gh.cells.columns(ph(c_fine)); % column index for fine cells part of top surface subgrid
      ve_col = Gh.cells.columns(ph(c_ve)); % column index for ve cells part of top surface subgrid
      top_surface_col = Gh.cells.columns(ph(Gt.cells.global)); % column index of cells definint top surface
      
      %ctrap_fine = ismember(Gt.cells.global, c_fine); % .global returns global index from full-dim grid
      ctrap_fine = ismember(top_surface_col, fine_col);
      ctrap_ve = ismember(top_surface_col, ve_col);
    
      % NB: only select top surface fine cells for ve regions. For fine
      % regions we need ALL fine cells
      kkf = kk(c_ve);                 
      kkf = c_ve(kkf == min(kkf));          
      cf_ve = c_ve(ismember(c_ve, kkf)); % only select topmost fine cells of associated ve columns from hybrid grid       
      ch_ve = ph(cf_ve); % should be same size as traps_ve, since only top (fine) cells are chosen
      ch_fine = ph(c_fine); % parent selection for ALL cells in fine region of subgrid
    
    [masses_ve, masses_ve0] = massTrappingDistributionVE(Gt, Gh, ch_ve, p, sG, sW, h, h_max, ...
                                                                 h_B, h_H, h_BH, Hh, Hbh, cB, cH, cBH, ...
                                                                 rock, fluidADI, traps, ctrap_ve, dh, opt, varargin{:});
    % The subgrid under the top surface might include fine cells -> these
    % cells require separate trap analysis
    [masses_fine, masses_fine0] = massTrappingDistributionFine(Gt, Gh, ch_fine, p, sG, sW, ...
                                                                 rock, fluidADI, traps, ctrap_fine, dh, opt, varargin{:});
                                                             
    % sum up masses for ve parts and fine parts of current subgrid
     masses = masses_ve + masses_fine;
     masses0 = cell2mat(masses_ve0) + cell2mat(masses_fine0);   
end


function [masses, masses_0] = massTrappingDistributionVE(Gt, Gh, c_ve, p, sG, sW, h, h_max, ...
                                                                 h_B, h_H, h_BH, Hh, Hbh, cB, cH, cBH, ...
                                                                 rock, fluidADI, traps, ctraps_ve, dh, opt, varargin)
    % DESCR: Compute masses for VE regions bounded by a sealing layer
    % Need to account for residual content originating from veBottom and
    % veHorizontal transition regions.
    ct_hybrid_height = Gh.cells.height(c_ve);
    %ct_hybrid_z = (Gh.cells.topDepth(c_ve) + Gh.cells.bottomDepth(c_ve))/2; % to get centroid
    ct_hybrid_z = Gh.cells.topDepth(c_ve);
        
    [cB_idx, cB_sub, ...
          cH_idx, cH_sub, ...
          cBH_idx, cBH_sub] = getRelaxedVEIndices(c_ve, cB, cH, cBH);    
    
    % Extracting relevant information from 'sol'
    sw=fluidADI.krPts.w(1);%liquid residual saturation (scalar)
    sr=fluidADI.krPts.g(1);%gas residual saturation (scalar)
    SW = sW(c_ve) .* ct_hybrid_height; % used to check that all elements in trapping inventory sum up to net fluid amount
    SG = sG(c_ve) .* ct_hybrid_height; % height of mobile plume
    h = h(c_ve);
    h_max = h_max(c_ve);
    rs = opt.rs;
    p = p(c_ve);
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro(c_ve) .* Gt.cells.volumes(ctraps_ve) .* pvMult; % think Gt.cells.volumes is fine -> it's just the surface area of each top surface cell
    if isfield(rock,'ntg')
       pv = rock.poro(c_ve) .* Gt.cells.volumes(ctraps_ve) .* rock.ntg .* pvMult;
       %effective area accounting for possible net-to-gross data
    end
    rhoCO2   = fluidADI.rhoGS .* fluidADI.bG(p);
    gasPhase = sum(pv .* (rhoCO2 .* SG));
    
    if isfield(traps, 'z_spill_loc')
        % trap heights -> only calculated from geometry, not h or h_max
        zt = max(traps.z_spill_loc(ctraps_ve) - ct_hybrid_z, 0); % Gt.cells.z --> ct_hybrid_z
        zt = min(zt, ct_hybrid_height); % Gt.cells.H --> ct_hybrid_height
    else
        tr   = traps.trap_regions(ctraps_ve); 
        tr_u = setdiff(unique(tr), 0); % get unique trap regions (except spill locations, = 0) for VE parts of subgrid to query correct z-values for traps
        tr_z = [min(ct_hybrid_z); traps.trap_z(tr_u)]; % trap region zero -> spills out of domain -> set z-value of spill point to highest point in grid
        tr   = 1 + tr;
        tr(tr>numel(tr_z)) = 1; % to ensure all trap regions with no calculated heights are assigned zero spill region height
        zt   = max(tr_z(tr) - ct_hybrid_z, 0);
        zt   = min(zt, ct_hybrid_height);
    end
    
    % NB: subtraps not yet accounted for in hybrid framework!
    h_sub  = zeros(numel(c_ve), 1); % subtrapped part of 'h'
    hm_sub = zeros(numel(c_ve), 1); % subtrapped part of 'h_max'
    if ~isempty(dh)
        h_sub  = min(dh, h);
        hm_sub = min(dh, h_max);
    end
    h_eff  = h - h_sub;
    hm_eff = h_max - hm_sub;
    
    [h_res_struct, h_res_free] = computeTrappedHeightsRVE(zt, h_eff, hm_eff, ct_hybrid_height, ...
                                                           h_B, h_H, h_BH, ...
                                                           Hh, Hbh, ...
                                                           cB, cB_idx, cB_sub, ...
                                                           cH, cH_idx, cH_sub, ...
                                                           cBH, cBH_idx, cBH_sub);

    strucVol  = sum(min(zt, h_eff) .* pv .* rhoCO2);             % trapped, flowing
    plumeVol  = sum(rhoCO2 .* h_eff.* pv) - strucVol;            % non-trapped, flowing
    resStruc  = (strucVol + sum(h_res_struct .* rhoCO2 .* pv)) * sr;      % trapped, res
    freeStruc = strucVol * (1 - sr - sw);                          % trapped, non-res
    freeRes   = plumeVol * sr;                                     % non-trapped, flowing, res
    freeMov   = plumeVol * (1 - sw - sr);                          % non-trapped, flowing, non-res
    
    resTrap = sum(h_res_free.*rhoCO2.*pv) .* sr;
    %resTrap   = sum(max(hm_eff - max(zt, h_eff), 0) .* ...
    %                rhoCO2 .* pv ) .* sr;                          % non-trapped, non-flowing, res
    resDis    = fluidADI.rhoGS .* sum(pv.* (rs .* fluidADI.bW(p) .* SW)); % dissolved
    subtrap   = sum((hm_sub * sr + h_sub * (1 - sr - sw)) .* pv .* rhoCO2);
          
    masses    = max([value(resDis), value(resStruc), value(resTrap), ...
                     value(freeRes), value(freeStruc), value(subtrap), ...
                     value(freeMov)], 0); % may be ADI variables
         
    if(abs(sum(masses(2:end))-gasPhase) > 1e-3 * gasPhase)
        %disp('There is a mismatch between mass calculations');
        %fprintf('Mass: %d. Gas: %d\n', abs(sum(masses(2:end))-gasPhase), 1e-3 * gasPhase);       
    end
    
    if nargout > 1
        % values one per cell of grid Gt
        strucVol_0  = min(zt, h_eff) .* pv .* rhoCO2;
        plumeVol_0  = rhoCO2 .*h_eff .* pv - strucVol_0;
        resStruc_0  = (strucVol_0 + h_res_struct.*rhoCO2.*pv) .* sr;
        freeStruc_0 = strucVol_0 .* (1-sr-sw);
        freeRes_0   = plumeVol_0 .* sr;
        freeMov_0   = plumeVol_0 .* (1-sr-sw);
        resTrap_0   = max(hm_eff - max(zt, h_eff),0) .* rhoCO2 .* pv .* sr;
        resDis_0    = fluidADI.rhoGS .* pv .* (rs .* fluidADI.bW(p) .* SW);
        subtrap_0   = (hm_sub * sr + h_sub * (1 - sr - sw)) .* pv .* rhoCO2;

        masses_0 = [{resDis_0}, {resStruc_0}, {resTrap_0}, {freeRes_0}, ...
                    {freeStruc_0}, {subtrap_0}, {freeMov_0}]; % may be ADI vars
    end
                                                             
end

function [masses, masses_0] = massTrappingDistributionFine(Gt, Gh, c_fine, p, sG, sW, ...
                                                             rock, fluidADI, traps, ctraps_fine, dh, opt, varargin)
    opt.rs = 0;
    opt    = merge_options(opt, varargin{:});        
    
    ct_hybrid_height = Gh.cells.height(c_fine);
    ct_hybrid_top = Gh.cells.topDepth(c_fine);
    ct_hybrid_bottom = Gh.cells.bottomDepth(c_fine);
    ct_hybrid_centroid = (ct_hybrid_top + ct_hybrid_bottom)/2;
    cols = Gh.cells.columns(c_fine);
    %[col, idx] = sort(cols);
    cols_u = unique(cols);
        
    top_z = accumarray(cols, ct_hybrid_centroid, [], @min); % global top of fine region
    bottom_z = accumarray(cols, ct_hybrid_centroid, [], @max); % global bottom of fine region
    top_z = top_z(cols_u);
    bottom_z = bottom_z(cols_u);
    tot_height = bottom_z - top_z;
    
    % Extracting relevant information from 'sol'
    sw=fluidADI.krPts.w(1);%liquid residual saturation (scalar)
    sr=fluidADI.krPts.g(1);%gas residual saturation (scalar)
    SW = sW(c_fine); % used to check that all elements in trapping inventory sum up to net fluid amount
    SG = sG(c_fine); % height of mobile plume
    
    rs = opt.rs;
    p = p(c_fine);
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro(c_fine) .* Gh.cells.volumes(c_fine) .* pvMult; % since cells ae full-dim, use Gh.cells.volumes, not volumes from top surface
    if isfield(rock,'ntg')
       pv = rock.poro(c_fine) .* Gh.cells.volumes(c_fine) .* rock.ntg .* pvMult;
       %effective area accounting for possible net-to-gross data
    end
    rhoCO2   = fluidADI.rhoGS .* fluidADI.bG(p);
    gasPhase = sum(pv .* (rhoCO2 .* SG));        
    
    
    if isfield(traps, 'z_spill_loc')
        % trap heights -> only calculated from geometry, not h or h_max
        zt = traps.z_spill_loc(ctraps_fine); % this is a DEPTH, not a THICKNESS
    else
        tr   = traps.trap_regions(ctraps_fine);
        tr_z = [min(top_z); traps.trap_z];
        tr   = 1 + tr;     
        zt = tr_z(tr); % DEPTH, not THICKNESS
    end
    
    z = Gh.cells.centroids(:,3);
    z = z(c_fine);  
    
    freeStruc = 0; resStruc = 0; freeRes = 0;
    freeMov = 0; resTrap = 0; resDis = 0;
    
    for i=1:numel(cols_u)
        cols_i = cols_u(i);
        col_idx = cols == cols_i; % all cells connected to cols_i
        zi = z(col_idx); % extract depths for all fine cells in this ve column
        SG_i = SG(col_idx);
        SW_i = SW(col_idx);
        p_i = p(col_idx);
        pv_i = pv(col_idx);
        rho_pv = fluidADI.rhoGS.*pv_i;
                     
        strucTrapped_i = zi < zt(i);
        % --- RESIDUAL BUFFER --- (to separate residual and mobile CO2)
        res_buff = 1.5; 
        % -----------------------
        freePlume_i = zi >= zt(i) & SG_i > res_buff*sr; % cells in mobile plume
        resPlume_i = zi >= zt(i) & SG_i <= res_buff*sr; % cells in immobilized plume
         
        freeStruc = freeStruc + sum((max(SG_i, sr) - sr).*rho_pv.* strucTrapped_i); % sum over all fine cells in column that ae structurally trapped
        resStruc = resStruc + sum(min(SG_i, sr).*rho_pv .* strucTrapped_i);
        freeRes = freeRes + sum(sr.*rho_pv .* freePlume_i); % we know SG_i > sr, so residual value can be fixed at sr
        freeMov = freeMov + sum((SG_i - sr).*rho_pv .* freePlume_i);
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

        masses_0 = [{resDis_0}, {resStruc_0}, {resTrap_0}, {freeRes_0}, ...
                    {freeStruc_0}, {subtrap_0}, {freeMov_0}]; % may be ADI vars
    end
end


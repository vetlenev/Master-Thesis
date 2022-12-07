function [masses, masses_0] = massTrappingFine_TopSurface(Gh, cmapf, p, sG, sW, ...
                                                            rock, fluidADI, trapstruct, varargin)
% Compute the trapping status distribution of CO2 in fine regions of global
% top-surface grid, removing cells part of top surface of semi-perm layers
% NOTE: Top-surface grid Gt not needed as argument because volumes of fine
% cells can be directly extracted from hybrid grid Gh.
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
%   ?Gti?        - cell array of semi-permeable top surface grids
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
    
    ct_hybrid_height = Gh.cells.height(cmapf);
    ct_hybrid_top = Gh.cells.topDepth(cmapf);
    ct_hybrid_bottom = Gh.cells.bottomDepth(cmapf);
    cols = Gh.cells.columns(cmapf);
    %[col, idx] = sort(cols);
    cols_u = unique(cols);
        
    top_z = accumarray(cols, ct_hybrid_top, [], @min); % global top of fine region
    bottom_z = accumarray(cols, ct_hybrid_bottom, [], @max); % global bottom of fine region
    top_z = top_z(cols_u);
    bottom_z = bottom_z(cols_u);
    tot_height = bottom_z - top_z;
    
    % Extracting relevant information from 'sol'
    sw=fluidADI.krPts.w(1);%liquid residual saturation (scalar)
    sr=fluidADI.krPts.g(1);%gas residual saturation (scalar)
    SW = sW(cmapf); % used to check that all elements in trapping inventory sum up to net fluid amount
    SG = sG(cmapf); % height of mobile plume
    
    rs = opt.rs;
    p = p(cmapf);
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro(cmapf) .* Gh.cells.volumes(cmapf) .* pvMult; % think Gt.cells.volumes is fine -> it's just the surface area of each top surface cell
    if isfield(rock,'ntg')
       pv = rock.poro(cmapf) .* Gh.cells.volumes(cmapf) .* rock.ntg .* pvMult;
       %effective area accounting for possible net-to-gross data
    end
    rhoCO2   = fluidADI.rhoGS .* fluidADI.bG(p);
    gasPhase = sum(pv .* (rhoCO2 .* SG));
    
    
    
    
    if isfield(trapstruct, 'z_spill_loc')
        % trap heights -> only calculated from geometry, not h or h_max
        zt = max(trapstruct.z_spill_loc - top_z, 0); % Gt.cells.z --> ct_hybrid_z
        zt = min(zt, tot_height); % Gt.cells.H --> ct_hybrid_height
    else
        tr   = trapstruct.trap_regions; 
        tr_z = [min(top_z); trapstruct.trap_z]; % trap region zero -> spills out of domain -> set z-value of spill point to highest point in grid
        tr   = 1 + tr;
        tr(tr>numel(tr_z)) = 1; % to ensure all trap regions with no calculated heights are assigned zero spill region height
        zt   = max(tr_z(tr) - top_z, 0);
        zt   = min(zt, tot_height);
    end
    
    z = Gh.cells.centroids(:,3);
    z = z(cmapf);  
    
    freeStruc = 0; resStruc = 0; freeRes = 0;
    freeMov = 0; resTrap = 0; resDis = 0;
    
    for i=1:numel(cols_u)
        cols_i = cols_u(i);
        col_idx = cols == cols_i;
        zi = z(col_idx);
        SG_i = SG(col_idx);
        SW_i = SW(col_idx);
        p_i = p(col_idx);
        pv_i = pv(col_idx);
        rho_pv = fluidADI.rhoGS.*pv_i;
        % Compare: zi < zt(i) ?
        strucTrapped_i = zi < zt(i);
        freePlume_i = zi >= zt(i) & SG_i > sr;
        resPlume_i = zi >= zt(i) & SG_i <= sr;
        
        freeStruc = freeStruc + sum((max(SG_i, sr) - sr).*rho_pv.* strucTrapped_i);
        resStruc = resStruc + sum(min(SG_i, sr).*rho_pv .* strucTrapped_i);
        freeRes = freeRes + sum(sr.*rho_pv .* freePlume_i); % we know SG_i > sr, so residual value can be fixed at sr
        freeMov = freeMov + sum((max(SG_i, sr) - sr).*rho_pv .* freePlume_i);
        resTrap = resTrap + sum(min(SG_i, sr).*rho_pv .* freePlume_i);
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
        strucVol_0  = min(zt, h_eff) .* pv .* rhoCO2;
        plumeVol_0  = rhoCO2 .*h_eff .* pv - strucVol_0;
        resStruc_0  = (strucVol_0 + hdift.*rhoCO2.*pv) .* sr;
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

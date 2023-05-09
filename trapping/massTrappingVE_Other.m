function [masses, masses_0] = massTrappingVE_Other(Gh, cVE_h, p, sG, sW, h, h_max, ...
                                                     h_B, h_H, h_BH, Hh, Hbh, cB, cH, cBH, ...
                                                     rock, fluidADI, varargin)
% Compute the trapping status distribution of CO2 in each cell of global
% top-surface grid, removing cells part of top surface of semi-perm layers
%
% DESCRIPTION:
%
% PARAMETERS:
%   Gh         - hybrid grid
%   cVE_h      - desired cells with indices from hybrid grid
%   p          - pressure, one value per cell of grid
%   sW         - water saturation, one value per cell of HYBRID grid
%   sG         - gas saturation, one value per cell of grid
%   h          - gas height below top surface, one value per cell of grid
%   h_max      - maximum historical gas height
%   h_B        - height of residual CO2 originating from bottom sealing
%                layer
%   h_H        - height of residual CO2 originating from horizontal sealing
%                layer(s)
%   h_BH       - height of residual CO2 originating from bottom AND horizontal sealing
%                layers
%   Hh         - height of horizontal sealing layers
%   Hbh        - height of bottom+horizontal sealing layers
%   cB         - bottom sealing cells
%   cH         - horizontal sealing cells
%   cBH        - bottom+horizontal sealing cells
%   rock       - rock parameters corresponding to 'Gt'
%   fluidADI   - ADI fluid object (used to get densities and compressibilities)
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

    cVE_h = unique(cVE_h);

    ct_hybrid_height = Gh.cells.height(cVE_h);

    [cB_idx, cB_sub, ...
      cH_idx, cH_sub, ...
      cBH_idx, cBH_sub] = getRelaxedVEIndices(cVE_h, cB, cH, cBH); 
    
    % Extracting relevant information from 'sol'
    sw=fluidADI.krPts.w(1);%liquid residual saturation (scalar)
    sr=fluidADI.krPts.g(1);%gas residual saturation (scalar)
    SW = sW(cVE_h) .* ct_hybrid_height;
    SG = sG(cVE_h) .* ct_hybrid_height;
    h = h(cVE_h);
    h_max = h_max(cVE_h);
    rs = opt.rs;    
    p = p(cVE_h); 
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro(cVE_h) .* Gh.cells.volumes(cVE_h)./ct_hybrid_height .* pvMult; % volumes ./ height to get lateral surface area of ve cell
    if isfield(rock,'ntg')
       pv = rock.poro(cVE_h) .* Gh.cells.volumes(cVE_h)./ct_hybrid_height .* rock.ntg(cVE_h) .* pvMult;
       %effective area accounting for possible net-to-gross data
    end
    rhoCO2   = fluidADI.rhoGS .* fluidADI.bG(p);                  
                      
    zt = zeros(numel(cVE_h),1); % min possible spill-point depth => no CO2 strucutrally trapped

    [~, h_res_free] = computeTrappedHeightsRVE(zt, h, h_max, ct_hybrid_height, ...
                                                           h_B, h_H, h_BH, ...
                                                           Hh, Hbh, ...
                                                           cB, cB_idx, cB_sub, ...
                                                           cH, cH_idx, cH_sub, ...
                                                           cBH, cBH_idx, cBH_sub);
    
    % No structural trapping
    resStruc  = 0;      % trapped, res
    freeStruc = 0;      % trapped, non-res
        
    plumeVol = sum(rhoCO2 .* h.* pv); % Entire h is part of plume now
    freeRes   = plumeVol * sr;
    freeMov   = plumeVol * (1 - sw - sr);
    resTrap = sum(h_res_free.*rhoCO2.*pv) .* sr;
    
    resDis    = fluidADI.rhoGS .* sum(pv.* (rs .* fluidADI.bW(p) .* SW)); % dissolved
    subtrap   = 0;
          
    masses    = max([value(resDis), value(resStruc), value(resTrap), ...
                     value(freeRes), value(freeStruc), value(subtrap), ...
                     value(freeMov)], 0); % may be ADI variables            
    
    if nargout > 1
        % values one per cell of grid Gt                    
        resStruc_0  = zeros(numel(SG), 1);
        freeStruc_0 = zeros(numel(SG), 1);

        plumeVol_0 = rhoCO2 .* h.* pv; % Entire h is part of plume now
        freeRes_0   = plumeVol * sr;
        freeMov_0   = plumeVol * (1 - sw - sr);
        resTrap_0 = h_res_free.*rhoCO2.*pv .* sr;
        
        resDis_0    = fluidADI.rhoGS .* pv .* (rs .* fluidADI.bW(p) .* SW);
        subtrap_0   = 0;

        masses_0 = [{resDis_0}, {resStruc_0}, {resTrap_0}, {freeRes_0}, ...
                    {freeStruc_0}, {subtrap_0}, {freeMov_0}]; % may be ADI vars
    end
end

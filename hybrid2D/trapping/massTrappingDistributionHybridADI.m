function [masses, masses_0] = massTrappingDistributionHybridADI(Gt, Gh, cTopHybrid, p, sG, sW, h, h_max, ...
                                                                 h_B, h_H, h_BH, Hh, Hbh, cB, cH, cBH, ...
                                                                 rock, fluidADI, trapstruct, dh, varargin)
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
    
    % cTopHybrid: fine cells at top surface of grid Gt with hybrid indices    
    cth = cTopHybrid;
    
    ct_hybrid_height = Gh.cells.height(cth);
    ct_hybrid_z = Gh.cells.topDepth(cth);
        
    cB_idx = ismember(cB, cth); % index for arrays with one element per bottom cell cB (from hybrid indexes for full grid)
    cB_subidx = ismember(cth, cB); % index for arrays with one element per cell in subgrid (from hybrid indexes for full grid)
    %cB = cB(cB_idx); % only select bottom VE flux cells that are part of subgrid Gt
    cB_sub = find(cB_subidx);
           
    cH_idx = ismember(cH, cth); % with duplicates    
    cH_subidx = ismember(cth, cH);
    cH = cH(cH_idx);
    cH_sub = find(cH_subidx); % indices for cHorz cells in subgrid
    [cH_u, ~, uHidx] = unique(cH); 
    
    if ~isempty(cH)
        cH_sub = extendIndexForDuplicates(cH, cH_u, cH_sub);
    end
        
    cBH_idx = ismember(cBH, cth);
    cBH_subidx = ismember(cth, cBH);
    cBH = cBH(cBH_idx); 
    cBH_sub = find(cBH_subidx);
    [cBH_u, ~, uBHidx] = unique(cBH);
    
    if ~isempty(cBH)
        cBH_sub = extendIndexForDuplicates(cBH, cBH_u, cBH_sub);
    end
    
    % Extracting relevant information from 'sol'
    sw=fluidADI.krPts.w(1);%liquid residual saturation (scalar)
    sr=fluidADI.krPts.g(1);%gas residual saturation (scalar)
    SW = sW(cth) .* ct_hybrid_height; % used to check that all elements in trapping inventory sum up to net fluid amount
    SG = sG(cth) .* ct_hybrid_height; % height of mobile plume
    h = h(cth);
    h_max = h_max(cth);
    rs = opt.rs;
    p = p(cth);
    
    pvMult = 1; 
    if isfield(fluidADI, 'pvMultR')
        pvMult =  fluidADI.pvMultR(p);
    end
    pv       = rock.poro(cth) .* Gt.cells.volumes .* pvMult; % think Gt.cells.volumes is fine -> it's just the surface area of each top surface cell
    if isfield(rock,'ntg')
       pv = rock.poro(cth) .* Gt.cells.volumes .* rock.ntg .* pvMult;
       %effective area accounting for possible net-to-gross data
    end
    rhoCO2   = fluidADI.rhoGS .* fluidADI.bG(p);
    gasPhase = sum(pv .* (rhoCO2 .* SG));
    
    if isfield(trapstruct, 'z_spill_loc')
        % trap heights -> only calculated from geometry, not h or h_max
        zt = max(trapstruct.z_spill_loc - ct_hybrid_z, 0); % Gt.cells.z --> ct_hybrid_z
        zt = min(zt, ct_hybrid_height); % Gt.cells.H --> ct_hybrid_height
    else
        tr   = trapstruct.trap_regions; 
        tr_z = [min(ct_hybrid_z); trapstruct.trap_z]; % trap region zero -> spills out of domain -> set z-value of spill point to highest point in grid
        tr   = 1 + tr;
        tr(tr>numel(tr_z)) = 1; % to ensure all trap regions with no calculated heights are assigned zero spill region height
        zt   = max(tr_z(tr) - ct_hybrid_z, 0);
        zt   = min(zt, ct_hybrid_height);
    end
    
    % NB: subtraps not yet accounted for in hybrid framework!
    h_sub  = zeros(Gt.cells.num, 1); % subtrapped part of 'h'
    hm_sub = zeros(Gt.cells.num, 1); % subtrapped part of 'h_max'
    if ~isempty(dh)
        h_sub  = min(dh, h);
        hm_sub = min(dh, h_max);
    end
    h_eff  = h - h_sub;
    hm_eff = h_max - hm_sub;
    
    % --- Heights of remaining residual plumes --- 
    % bottom: (h_B_all has one element per per VE cell in subgrid => index
    % by cB_sub. h_B has one element per bottom cells cB of full grid =>
    % index by cB_idx
    h_B_struct = zeros(size(zt));
    h_B_free = zeros(size(zt));
    % hB = ct_hybrid_height (i.e. height of VE cell in subgrid)   
    h_B_struct(cB_sub) = zt(cB_sub) - min(zt(cB_sub), h_B(cB_idx)) - max(zt(cB_sub) - ct_hybrid_height(cB_sub), 0);   
    h_B_free(cB_sub) = max(ct_hybrid_height(cB_sub), zt(cB_sub)) - max(h_B(cB_idx), zt(cB_sub));
    
    % horizontal:
    h_H_struct = zeros(size(zt)); % heights of residual plume (for cHorz) inside structural traps
    h_H_free = zeros(size(zt)); % heights of residual plume (for cHorz) outside structural traps
    h_H_structi = zt(cH_sub) - min(zt(cH_sub), h_H(cH_idx)) - max(zt(cH_sub) - Hh(cH_idx), 0); % heights for individual semi-perm layers
    h_H_freei = max(Hh(cH_idx), zt(cH_sub)) - max(h_H(cH_idx), zt(cH_sub));

    % bottom + horizontal:
    h_BH_struct = zeros(size(zt)); % one fo each cell in subgrid - multiple values in cHorz cells to be accumulated
    h_BH_free = zeros(size(zt)); % heights of residual plume (for cHorz+cBottom) outside structural traps
    h_BH_structi = zt(cBH_sub) - min(zt(cBH_sub), h_BH(cBH_idx)) - max(zt(cBH_sub) - Hbh(cBH_idx), 0);        
    h_BH_freei = max(Hbh(cBH_idx), zt(cBH_sub)) - max(h_BH(cBH_idx), zt(cBH_sub));

    %h_res_free = max(hm_eff - max(h, zt), 0); % residual heights outside structural traps
    % sum up residual plume heights inside structural traps    
    if ~isempty(cH)
        [cH_sub_u,~,cH_sub_map] = unique(cH_sub);
        
        h_H_sum.struct = accumarray(cH_sub_map, h_H_structi); % use c to accumulate for correct index
        h_H_struct(cH_sub_u) = h_H_sum.struct;   
        
        h_H_sum.free = accumarray(cH_sub_map, h_H_freei);       
        h_H_free(cH_sub_u) = h_H_sum.free; % add height of additional residual plumes       
    end
    
    if ~isempty(cBH)
        [cBH_sub_u,~,cBH_sub_map] = unique(cBH_sub);
        
        h_BH_sum.struct = accumarray(cBH_sub_map, h_BH_structi); % use c to accumulate for correct index
        h_BH_struct(cBH_sub_u) = h_BH_sum.struct;  
        
        h_BH_sum.free = accumarray(cBH_sub_map, h_BH_freei);
        h_BH_free(cBH_sub_u) = h_BH_sum.free;
    end
    % ------------------------------------------
    h_res_free = max(hm_eff - max(h, zt), 0) + ... % top residual plume below trap
                 h_B_free + h_H_free + h_BH_free; % remaining residual plumes below trap
    
    % this requires that the fluid has a sharp interface relperm of normal type    
    hdift     = max(min(zt, hm_eff) - min(zt, h_eff),0);    % trapped part of h_max-h INSIDE structural trap
    % --- Handle remaining residual plumes ---
    % Not affected by subtraps ??  
    hdift_B = max(h_B_struct, 0);          
    hdift_H = max(h_H_struct, 0);      
    hdift_BH = max(h_BH_struct, 0);
    hdift_net = hdift + hdift_B + hdift_H + hdift_BH;
    
    strucVol  = sum(min(zt, h_eff) .* pv .* rhoCO2);             % trapped, flowing
    plumeVol  = sum(rhoCO2 .* h_eff.* pv) - strucVol;            % non-trapped, flowing
    resStruc  = (strucVol + sum(hdift_net .* rhoCO2 .* pv)) * sr;      % trapped, res
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

function [cH_sub] = extendIndexForDuplicates(cH, cH_u, cH_sub)
    cH_ctf = [cH_u, cH_sub]; % make matrix mapping from hybrid cell indices to associated local index of subgrid    
    cH_ctf1 = cH_ctf(:,1);
    cH_ctf2 = cH_ctf(:,2);

    edges = min(cH):max(cH)+1; % +1 to include last index
    [counts, values] = histcounts(cH, edges);
    duplicates = values(counts > 1);
    counts = counts(counts > 1);        
    cnt = 0; % to count number of elements cH_subidx is extended with

    cH_sub = zeros(size(cH));
    idx_cH = ~ismember(cH, duplicates); % no duplicate indexes for cells cH
    cH_no_duplicates = cH(idx_cH);
    idx_ctf = ismember(cH_ctf1, cH_no_duplicates); % no duplicate indexes for cells ctf
    cH_sub(idx_cH) = cH_ctf2(idx_ctf); % add no-duplicate elements by vectorization
    
    for i=1:numel(duplicates) % extend cH_sub sequentially with repeated cells
        d = duplicates(i);
        %d_idx = find(ismember(cH_ctf1, d)) + cnt;
        d_idx_glob = ismember(cH, d); % logical index for duplicate
        d_sub = cH_ctf2(ismember(cH_ctf1, d)); % mapped value
        % CHECK IF THIS WORKS:        
        cH_sub(d_idx_glob) = d_sub; % insert mapped value into associated position in cH_sub
%         d_sub = cH_ctf2(ismember(cH_ctf1, d));
%         cnt = cnt + counts(i)-1;
%         cH_sub = [cH_sub(1:d_idx); repmat(d_sub, 1, counts(i)-1); cH_sub(d_idx+1:end)]; % add duplicates to subgrid indices
    end
end
function [h_res_struct, h_res_free] = computeTrappedHeightsRVE(zt, h, hmax, rve_height, ...
                                                           h_B, h_H, h_BH, ...
                                                           Hh, Hbh, ...
                                                           cB, cB_idx, cB_sub, ...
                                                           cH, cH_idx, cH_sub, ...
                                                           cBH, cBH_idx, cBH_sub)
% Compute heights of remaining residual plumes, i.e., relaxed VE columns,
% separating parts that are structurally trapped and part of the free
% plume.
% PARAMTERS:
%   zt      - spill point depth of the trap a cell is associated with
%   h       - height of top mobile plume
%   hmax    - max reached height of top mobile plume
%   rve_height - height of RVE column
%   h_X     - depth from top of VE column to top of plume rising from
%               relaxed VE cell X.
%   Hx      - height of virtual cell associated with relaxed VE cell x.
%   cX      - all hybrid indices of relaxed VE column of type X
%   cX_idx  - hybrid indices of relaxed VE columns of type X from a subset
%               of cells from full grid 
%   cX_sub  - (local) subgrid indices of relaxed VE columns of type X.
% RETURNS:
%   h_res_struct - net height of residual CO2 structurally trapped
%   h_res_free - net height of residual CO2 below traps

    h_B_struct = zeros(size(zt));
    h_B_free = zeros(size(zt));
    % hB = ct_hybrid_height (i.e. height of VE cell in subgrid)   
    h_B_struct(cB_sub) = zt(cB_sub) - min(zt(cB_sub), h_B(cB_idx)) - max(zt(cB_sub) - rve_height(cB_sub), 0);   
    h_B_free(cB_sub) = max(rve_height(cB_sub), zt(cB_sub)) - max(h_B(cB_idx), zt(cB_sub));
    
    % horizontal:
    h_H_struct = zeros(size(zt)); % heights of residual plume (for cHorz) inside structural traps
    h_H_free = zeros(size(zt)); % heights of residual plume (for cHorz) outside structural traps
    h_H_structi = zt(cH_sub) - min(zt(cH_sub), h_H(cH_idx)) - max(zt(cH_sub) - Hh(cH_idx), 0); % heights for individual semi-perm layers
    h_H_freei = max(Hh(cH_idx), zt(cH_sub)) - max(h_H(cH_idx), zt(cH_sub));

    % bottom + horizontal:
    h_BH_struct = zeros(size(zt)); % one for each cell in hybrid grid - multiple values in cHorz cells to be accumulated
    h_BH_free = zeros(size(zt)); % heights of residual plume (for cHorz+cBottom) outside structural traps
    h_BH_structi = zt(cBH_sub) - min(zt(cBH_sub), h_BH(cBH_idx)) - max(zt(cBH_sub) - Hbh(cBH_idx), 0);        
    h_BH_freei = max(Hbh(cBH_idx), zt(cBH_sub)) - max(h_BH(cBH_idx), zt(cBH_sub));
   
    % sum up residual plume heights inside structural traps    
    if ~isempty(cH_sub)
        [cH_sub_u,~,cH_sub_map] = unique(cH_sub);
        
        h_H_sum.struct = accumarray(cH_sub_map, h_H_structi); % use duplicated index cH_sub_map to accumulate for correct index
        h_H_struct(cH_sub_u) = h_H_sum.struct;   
        
        h_H_sum.free = accumarray(cH_sub_map, h_H_freei);       
        h_H_free(cH_sub_u) = h_H_sum.free; % add height of additional residual plumes       
    end
    
    if ~isempty(cBH_sub)
        [cBH_sub_u,~,cBH_sub_map] = unique(cBH_sub);
        
        h_BH_sum.struct = accumarray(cBH_sub_map, h_BH_structi);
        h_BH_struct(cBH_sub_u) = h_BH_sum.struct;  
        
        h_BH_sum.free = accumarray(cBH_sub_map, h_BH_freei);
        h_BH_free(cBH_sub_u) = h_BH_sum.free;
    end
    % ------------------------------------------
    h_res_free = max(hmax - max(h, zt), 0) + ... % top residual plume below trap
                 h_B_free + h_H_free + h_BH_free; % remaining residual plumes below trap
    
    % this requires that the fluid has a sharp interface relperm of normal type    
    h_struct_res_top     = max(min(zt, hmax) - min(zt, h),0);    % trapped part of h_max-h (top plume) INSIDE structural trap
    % --- Handle remaining residual plumes ---
    % Not affected by subtraps ??  
    h_struct_B = max(h_B_struct, 0);          
    h_struct_H = max(h_H_struct, 0);      
    h_struct_BH = max(h_BH_struct, 0);
    h_res_struct = h_struct_res_top + h_struct_B + h_struct_H + h_struct_BH; % net trapped CO2 volumes inside structural traps -> accounts for parts of every plume that resides within a trap    
end
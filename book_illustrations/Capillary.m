classdef Capillary
    %CAPILLARY Functions for capillary pressure
    
    methods (Static)
        function pc_val = PcGas(Sg, swr, snr, p_e, cap, n)
            Sw = 1 - Sg;
            S_scaled = max( (Sw-swr)/(1-snr-swr), 1e-5);
            pc_val = p_e*S_scaled.^(-1/n); % Corey model
            pc_val(pc_val>cap) = cap; % Cap to prevent infinity
            pc_val(Sw<=swr) = cap;
            pc_val(Sw==1) = 0; % No pressure if no saturation
        end  
      
      function pc_val = PcNew(S, swr, snr, p_e, cap, n)
         S_scaled = max( (S-swr)/(1-snr-swr), 1e-5);
         pc_val = p_e*S_scaled.^(-1/n); % Corey model
         pc_val(pc_val>cap) = cap; % Cap to prevent infinity
         pc_val(S<=swr) = cap;
         pc_val(S==1) = 0; % No pressure if no saturation
      end  
      
      function pc_val = LeverettJ(S, phi, K, K_base, median_pc)                 
         %S_scaled = max((S - swr) / (1 - snr- swr), 1e-5); % NB: no scaling of saturation for Leverett-J                                   
         surf_tension = median_pc/sqrt(median(phi)/median(K)); % median reservoir properties give cap pressure of median_pc (barsa)
         pc_val = surf_tension*sqrt(phi./K).*min(max((1 - S), 0), 1);         
         pc_val(K == K_base) = 0; % No capillary pressure in background
      end
      
  
      function pc = runStandardPc(S, dummy_S, swr, snr, p_e, p_cap, layers, G)
        pc_vals = Capillary.PcGas(dummy_S, swr, snr, p_e, p_cap, 2);

        region_table = {[dummy_S, zeros(numel(dummy_S), 1)], [dummy_S, pc_vals]}; % container for pc values in each region      
        region_idx = {setdiff(G.cells.indexMap, layers), layers};      
        
        pc = interpReg(region_table, S, region_idx);
      end
      
      function pc = runHybridPc(S, dummy_S, swr, snr, p_e, p_cap, layers, G, nn_fine)
        pc_vals = Capillary.PcGas(dummy_S, swr, snr, p_e, p_cap, 2);
        region_table = {[dummy_S, zeros(numel(dummy_S), 1)], [dummy_S, pc_vals]}; % container for pc values in each region
    
        if numel(nn_fine) == G.cells.num % evaluate for entire region     
            if isfield(G, 'partition') % fetch coarse representation for each fine cell
                coarseIndexMap = (1:G.cells.num)';
                region_idx = {setdiff(coarseIndexMap, layers), layers}; % region to interpolate (rest, lowperm)
            else
                region_idx = {setdiff(G.cells.indexMap, layers), layers};
            end
            %pc = interpReg(region_table, S, region_idx);
        else % consider specific transition regions            
            % NEED TO CONSIDER INDICES FOR ENTIRE DOMAIN, OR EXTRACT CORRESPONDING INDICES FOR S !!!
            % -----
            %S(nn_fine) or something ?
            % -----
            % Or take cells c (which S is extracted from) as input argument ??
            sealingIdx = ismember(nn_fine, layers);
            sealing = nn_fine(sealingIdx);
            noSealing = nn_fine(~sealingIdx);            
            region_idx = {noSealing, sealing};
            
            % Alternative B            
        end
        pc = interpReg(region_table, S, region_idx);
      end
      
      function pc = runLeverettJ_2D(S, dummy_S, phi, K, dummy_K, K_base, layers, G)
        [grid_Sw, grid_K] = ndgrid(dummy_S, dummy_K);
        pc_vals = Capillary.LeverettJ(grid_Sw, phi, grid_K, K_base);

        region_table = {{grid_Sw, grid_K, zeros(size(grid_Sw))}, ...
                         {grid_Sw, grid_K,  pc_vals}}; % container for pc values in each region
        region_idx = {setdiff(G.cells.indexMap, layers).', layers}; % region to interpolate (rest, lowperm)
        pc = interpReg2D(region_table, S, K, region_idx);
      end
    end
end


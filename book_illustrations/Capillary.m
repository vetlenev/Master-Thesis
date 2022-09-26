classdef Capillary
    %CAPILLARY Functions for capillary pressure
    
    methods (Static)
        function pc_val = PcGas(Sg, swr, snr, p_e, cap, n)
            Sw = 1 - Sg;
            S_scaled = max( (Sw-swr)/(1-snr-swr), 1e-5);
            pc_val = p_e*S_scaled.^(-1/n); % Corey model
            pc_val(pc_val>cap) = cap; % Cap to prevent infinity
            pc_val(Sw<swr) = cap;
            pc_val(value(Sw)==1) = 0; % No pressure if no saturation (cast to double, == not supported for ADI)
        end  
      
      function pc_val = PcNew(S, swr, snr, p_e, cap, n)
         S_scaled = max( (S-swr)/(1-snr-swr), 1e-5);
         pc_val = p_e*S_scaled.^(-1/n); % Corey model
         pc_val(pc_val>cap) = cap; % Cap to prevent infinity
         pc_val(S<swr) = cap;
         pc_val(S==1) = 0; % No pressure if no saturation
      end  
      
      function pc_val = LeverettJ(S, phi, K, K_base, median_pc)                 
         %S_scaled = max((S - swr) / (1 - snr- swr), 1e-5); % NB: no scaling of saturation for Leverett-J                                   
         surf_tension = median_pc/sqrt(median(phi)/median(K)); % median reservoir properties give cap pressure of median_pc (barsa)
         pc_val = surf_tension*sqrt(phi./K).*min(max((1 - S), 0), 1);         
         pc_val(K == K_base) = 0; % No capillary pressure in background
      end
      
  
      function pc = runStandardPc(S, dummy_S, swr, snr, pe_sealing, pe_rest, p_cap, layers, G)
        pc_sealing = Capillary.PcGas(dummy_S, swr, snr, pe_sealing, p_cap, 4);
        pc_rest = Capillary.PcGas(dummy_S, swr, snr, pe_sealing, p_cap, 4);
        %pc_rest = zeros(numel(dummy_S), 1);
        
        region_table = {[dummy_S, pc_rest], [dummy_S, pc_sealing]}; % container for pc values in each region      
        region_idx = {setdiff(G.cells.indexMap, layers), layers};      
        
        pc = interpReg(region_table, S, region_idx);
      end
          
      
      function pc = runHybridPc(S, swr, snr, pe_sealing, pe_rest, p_cap, n_sealing, isVE)           
        % For (fine) sealing cells, use Brooks-Corey with high entry pressure
        pc = Capillary.PcGas(S, swr, snr, pe_sealing, p_cap, 1.5);     
        % For VE columns, only include entry pressure at sharp interface
        % For fine cells (not sealing), use Brooks-Corey with lower entry pressure                
        if numel(isVE) > 1
            isVE = isVE(~n_sealing);    
        end        
        pc(~n_sealing) = isVE.*pe_rest + ~isVE.*Capillary.PcGas(S(~n_sealing), swr, snr, pe_rest, p_cap, 1.5); % First term: entry pressure at interface only nonzero if there actually is an interface        
        
      end % isVE.*pe_rest
            
    end
end


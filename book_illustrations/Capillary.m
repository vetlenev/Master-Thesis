classdef Capillary
    %CAPILLARY Functions for capillary pressure
    
    methods (Static)
        function pc_val = PcGas(Sg, swr, snr, p_e, n)
            Sw = 1 - Sg;
            cap = 5*p_e; % pc capped at 5 times entry pressure           
            S_scaled = @(Sw) max( (Sw-swr)/(1-snr-swr), 1e-5);
            pc_val = @(Sw) p_e*S_scaled(Sw).^(-1/n); % Corey model
            
            % Define lower extension to avoid infinity and ensure
            % continuity
            %pc_max = pc_val(SwMin);
            %pc_val(Sw <= SwMin) = pc_max .* exp(10.*(Sw-SwMin));
            
            pc_val = pc_val(Sw);
            pc_val(pc_val>cap) = cap; % Cap to prevent infinity
            pc_val(Sw<=swr) = cap;
            %pc_val(value(Sw) > 1-1e-9) = 0; % No pressure if no saturation (cast to double, == not supported for ADI)
        end  
      
      function pc_val = PcSharp(Sg, swr, snr, p_e, n)
         Sw = 1 - Sg;
         %cap = 5*p_e;
         cap = p_e;
         pc_val = repmat(p_e, numel(value(Sg)), 1); % Corey model
         %pc_val(pc_val>cap) = cap; % Cap to prevent infinity
         %pc_val(Sw<=swr) = cap;
         %pc_val(value(Sw) > 1-1e-9) = 0; % No pressure if no saturation
      end      
      
  
      function [pc] = runStandardPcSharp(S, dummy_S, swr, snr, pe_sealing, pe_rest, n_cells, G)
        % NB: CHANGE TO PcSharp IF RUNNING SHARP INTERFACE
        pc_sealing = Capillary.PcSharp(dummy_S, swr, snr, pe_sealing, 4);
        pc_rest = Capillary.PcSharp(dummy_S, swr, snr, pe_rest, 4);
        %pc_rest = zeros(numel(dummy_S), 1);
        
        region_table = {[dummy_S, pc_rest], [dummy_S, pc_sealing]}; % container for pc values in each region      
        if isfield(G.cells, 'indexMap') % fine grid           
            region_idx = {setdiff(G.cells.indexMap, n_cells), n_cells};
            pc = interpReg(region_table, S, region_idx);
        else % hybrid grid           
            n_sealing = n_cells(ismember(n_cells, G.sealingCells));
            idx_diff = setdiff(n_cells, n_sealing);
            not_sealing = ismember(n_cells, idx_diff);
            region_idx = {n_cells(not_sealing), n_sealing};   
            
            if numel(value(S)) == G.cells.num
                sG = S; % evaluated for entire domain
            else % only a subset of cells evaluated -> 
                sG = zeros(G.cells.num, 1);                    
                sG(region_idx{1}) = value(S(find(not_sealing))); % high-permeable regions
                sG(region_idx{2}) = value(S(find(~not_sealing))); % low-permeable regions
            end
            
            pc = interpReg(region_table, sG, region_idx);
            pc = pc(n_cells);
        end
               
      end
                
      function pc = runHybridPcSharp(S, swr, snr, pe_sealing, pe_rest, n_sealing, isVE)                  
        pc = Capillary.PcSharp(S, swr, snr, pe_sealing, 4);     
        % For VE columns, only include entry pressure at sharp interface
        % For fine cells (not sealing), use Brooks-Corey with lower entry pressure                     
        %isVE = isVE(~n_sealing);     
        %pc(~n_sealing) = isVE.*pe_rest.*(S(~n_sealing) >= 1e-10) ...
        %                + ~isVE.*Capillary.PcSharp(S(~n_sealing), swr, snr, pe_rest, 4);        
        pc(~n_sealing) = Capillary.PcSharp(S(~n_sealing), swr, snr, pe_rest, 4);
      end 
            
      function [pc] = SharpHysteresis(Sn, Sn_max, swr, snr, pe_sealing, pe_rest, n_sealing, isVE)
          pc = runHybridPcSharp(Sn, swr, 0, pe_sealing, pe_rest, n_sealing, isVE); % primary drainage
          pc(value(Sn_max) > 1e-7) = runHybridPcSharp(Sn, swr, snr, pe_sealing, pe_rest, n_sealing, isVE); % main drainage          
      end  
    end
end


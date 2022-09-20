classdef VolumeTrapping
    %VOLUMETRAPPING Functions for computing volume and trapping 
    
    methods (Static)
        function [tot_vol, region_vol] = Co2VolumeRatio(G, cells, S, rock, fluid)
          % Estimates volume of CO2 stored in provided cells.
          % Inputs:
          %     G: grid object
          %     cells: selected cells to estimate volume for (i.e. top)
          %     S: CO2 saturation values for each cell in grid at current state
          %     rock: rock object
          %     fluid: fluid object
          % Outputs:
          %     tot_vol: total volume of CO2 in grid
          %     region_vol: volume of CO2 in region represented by
          %     selected cells                   
          tot_vol = sum(G.cells.volumes.*S(G.cells.indexMap))*mean(rock.poro);%*fluid.rhoOS;               
          region_vol = sum(G.cells.volumes(cells).*S(cells))*mean(rock.poro);%*fluid.rhoOS;      
        end       
      
      function [res_vol] = Co2ResidualTrapped(G, S, sor, rock)
          % Estimates volume of CO2 residually trapped in domain.
          % Inputs:
          %     G: grid object         
          %     S: CO2 saturation values for each cell in grid at current state        
          %     sor: residual oil saturation        
          %     rock: rock object         
          % Outputs:
          %     res_vol: residual volume of CO2
          S_res = min(S, sor);         
          res_vol = sum(G.cells.volumes.*S_res)*mean(rock.poro);%*fluid.rhoOS;
      end
      
      
      function [struct_imperm, struct_lowperm, ...
                struct_res_imperm, struct_res_lowperm] = Co2StructuralTrapped(G, S, sor, trapped_imperm, trapped_lowperm, rock)
          % Estimates volume of CO2 structurally trapped in domain, defined
          % as all cells residing below concave shape but above bottommost
          % part of layer, including cells immediately below layer.
          % 
          % Inputs:
          %     G: grid object         
          %     S: CO2 saturation values for each cell in grid at current state
          %     swr: residual water saturation
          %     sor: residual oil saturation
          %     trapped_imperm: cell candidates definitely
          %     subject to structural trapping
          %     trapped_lowperm: cell candidates possibly subject to
          %     structral trapping, depending on if entry pressure reached
          %     rock: rock object         
          % Outputs:
          %     struct_imperm: permanently structurally trapped CO2 volume,
          %                     omitting residual component         
          %     struct_lowperm: temporary structurally trapped CO2 volume,
          %                         omitting residual component
          %     struct_res_imperm: residually trapped CO2 inside
          %                         imperm structural traps
          %     struct_res_lowperm: residually trapped CO2 inside lowperm
          %                         structural traps
                 
          non_trapped_cells = G.cells.indexMap(setdiff(G.cells.indexMap, horzcat(trapped_imperm, trapped_lowperm)));         
          S(non_trapped_cells) = 0; % omit cells not candidated for structural trapping
          
          S_imperm = S; S_lowperm = S;
          S_imperm(trapped_lowperm) = 0;
          S_lowperm(trapped_imperm) = 0;
          
          S_imperm = max(S_imperm-sor, 0); % omit residually trapped CO2
          S_lowperm = max(S_lowperm-sor, 0);
          struct_imperm = sum(G.cells.volumes.*S_imperm)*mean(rock.poro);                       
          struct_lowperm = sum(G.cells.volumes.*S_lowperm)*mean(rock.poro);
          
          % residually trapped in structural traps:
          S_imperm = sor.*(S_imperm >= sor);         
          S_lowperm = sor.*(S_lowperm >= sor);
          struct_res_imperm = sum(G.cells.volumes.*S_imperm)*mean(rock.poro);
          struct_res_lowperm = sum(G.cells.volumes.*S_lowperm)*mean(rock.poro);
      end
    end
end



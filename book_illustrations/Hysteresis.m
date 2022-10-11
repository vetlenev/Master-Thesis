classdef Hysteresis
    %HYSTERESIS To include hysteresis effects for CO2 drainage and water
    %imbibition
    
    methods (Static)                  
        
        function [krn] = Killough(Sn, Sni, Sni_max, Snr_max, krn_PD, krn_PI)             
            C = 1/Snr_max - 1/Sni_max; % assuming zero residual sat for primary drainage curve
            S_nr = Sni ./ (C*Sni + 1);          
            Sn_dot = Snr_max + ((Sn-S_nr).*(Sni_max-Snr_max))./(Sni-S_nr);
            % New imbibition curve
            krn = krn_PI(Sn_dot).*krn_PD(Sni)./krn_PD(Sni_max);  
            krn(Sni < 1e-3) = 0.0; % all cells with no CO2 history forced immobilized 
        end
                
    end
end


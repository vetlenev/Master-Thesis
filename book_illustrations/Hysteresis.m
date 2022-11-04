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
        
        function [pc] = KilloughPc(Sn, Sni, Sni_max, Snr_max, pc_PD, pc_PI)
            C = 1/Snr_max - 1/Sni_max;
            S_nr = Sni ./ (C*Sni + 1);          
            Sn_dot = Snr_max + ((Sn-S_nr).*(Sni_max-Snr_max))./(Sni-S_nr);
            % New imbibition curve
            pc = pc_PI(Sn_dot).*pc_PD(Sni)./pc_PD(Sni_max);  
            pc(Sni < 1e-3) = 0.0; % all cells with no CO2 history has not reached entry pressure 
        end                
        
        function [pc_scan] = Genuchten(Sw, SwMin, snr, swr, pc_func, sw_func, alpha, n, gamma_D, gamma_I)
            % Hysteretic capillary pressure based on van Genuchten model.
            % Includes primary and first order scanning curves (assumed to
            % be reversible so we end up at primary drainage after full
            % reversal).
            
            % INPUTS:
            %   Sw: water saturation to get capillary pressure for
            %   SwMin: minimum water sat (max gas sat) reached over history
            %       where SwMin >= swr
            %   snr: minimum residual gas sat
            %   swr: minimum residual water sat
            %   pc_func: specific Genuchten model for capillary pressure
            %   sw_func: inverse of pc_func, solved for sw
            %   alpha: inverse entry pressure parameter
            %   n: power factor
            %   gamma_D: scaling factor for drainage
            %   gamma_I: scaling factor for imbibition
            
            % RETURNS:
            %   pc: scanning curve emanating from turning point SwMin
            C = 1/snr - 1/(1-swr);
            sne = (1-SwMin)./(1+C.*(1-SwMin)); % s_gr_Delta in [niemi, 1998]
            
            pc = @(Sw, sne, gamma) pc_func(Sw, sne, swr, alpha, n, gamma); % van Genuchten pc curves
            
            pc_PD = @(Sw) pc(Sw, 0, gamma_D); % primary drainage curve
            pc_I = @(Sw) pc(Sw, sne, gamma_I); % imbibition curve emanating from sne
            
            Sw_PD = @(pc) sw_func(pc, 0, swr, alpha, n, gamma_D);
            Sw_PI = @(pc) sw_func(pc, snr, swr, alpha, n, gamma_I); % saturations for main imbibition curve
            Sw_I = @(pc) sw_func(pc, sne, swr, alpha, n, gamma_I); % saturations for wetting curve emanating from sne
            
            % pc at reversal Sw = SwMin (reversal happens at main drainage curve)
            pc_max = pc_PD(SwMin);
            SwMin_PD = SwMin; % saturation at primary drainage at reversal point is simply stored max saturation SnMax = 1-SwMin. Same as Sw_PD(pc_max);
            SwMin_I = Sw_I(pc_max); % use inverse sat func to get sat corresponding to primary imbibition pc value at reversal point
            % pc_PD(SwMin_PD) = pc_I(SwMin_I)
                     
            % get rescaled saturation argument at turning point sne
            Sw_scan = @(Sw, Sw_PD, Sw_I) Sw_I + (Sw - Sw_PD).*((1-sne-Sw_I)./(1-sne-Sw_PD));                                                        
            % get scanning curve emanating from reversal point using
            % main wetting curve with rescaled argument
            Sw_hat = Sw_scan(Sw, SwMin_PD, SwMin_I);
            pc_scan = pc_I(Sw_hat); % plug rescaled sat into expression for wetting curve emanating from sne          
        end
        
        function [pc] = pc_func(Sw, sne, swr, alpha, n, gamma)
            S_eff = (Sw-swr)./(1-swr-sne);
            pc = -1/(alpha^gamma) .*( S_eff.^(-1/((1 - 1/n)^gamma)) - 1 ).^(1/(n^gamma));
        end
        function [sw] = sw_func(pc, sne, swr, alpha, n, gamma)
            pc_scaled = (-pc.*(alpha.^gamma)).^(n^gamma);
            sw = swr + (1-swr-sne).*(pc_scaled + 1).^(-(1-1/n)^gamma);
        end
    end
end


function pW_local = getPwFromHeight(pt, t, b, h, h_max, pW, g, rhow, rhog, swr, snr, varargin)
    % Internal function: Reconstruct water pressure from VE assumption,
    % including RESDIUAL SATURATED ZONE.
    % T is top of column for each cell, t is top of cells to be
    % reconstructed, with b/B having the same interpretation for the
    % respective cell bottoms. h is the height of the gas plume in each
    % cell.
    % pW is the water pressure at the BOTTOM of the cell, g the gravity
    % magnitude and rhow/rhog the water and gas densities.

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
    opt = struct('ve_model', 'sharp_interface', 'cNVE', [], ...
                  'p_entry', 0, 'transition', 'vertical', ...
                  'res_type', []);
    opt = merge_options(opt, varargin{:});
    
    dz_below = min(max(b - (t+h_max), 0), b-pt); % water below pt
    h_below = max(t+h - pt, 0); % mobile CO2 below pt
    hmax_below = (b-pt) - h_below - dz_below; % residual CO2 below pt
       
    if isempty(opt.res_type)
        pW_local = pW - rhog.*g.*h_below - rhow.*g.*(hmax_below + dz_below);
    elseif strcmp(opt.res_type, 'full')
       pW_local = pW - rhog.*g.*((1-swr).*h_below + snr.*hmax_below) ...
                     - rhow.*g.*(swr.*h_below + (1-snr).*hmax_below + dz_below);
    elseif strcmp(opt.res_type, 'mobile')
        % Omit residual pressure, only mobile zones:
        pW_local = pW - rhog.*g.*(1-swr).*h_below ...
                      - rhow.*g.*((1-snr).*hmax_below + dz_below);
    end
              
    if any(opt.cNVE)
       % Special treatment for NVEs where residual sat filled from below
       cNVE = opt.cNVE;
             
       b_pt = b(cNVE) - pt(cNVE); % height from bottom to point pt
       hmax_below = max(min(h_max(cNVE), b_pt), 0); % max(x, 0) to ensure no negatives
       dz_below = b_pt - hmax_below; 
       
       if isempty(opt.res_type)
           pW_local(cNVE) = pW(cNVE) - rhow(cNVE).*g.*b_pt;
       elseif strcmp(opt.res_type, 'full')
            pW_local(cNVE) = pW(cNVE) - rhog(cNVE).*g.*snr.*hmax_below ...
                                      - rhow(cNVE).*g.*((1-snr).*hmax_below + dz_below);                                
       elseif strcmp(opt.res_type, 'mobile')
            % Omit residual pressure:
            pW_local(cNVE) = pW(cNVE) - rhow(cNVE).*g.*((1-snr).*hmax_below + dz_below); % only water is mobile in cNVE cells --> pressure contribution only from water
       end
    end
end
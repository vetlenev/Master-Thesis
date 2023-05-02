function pG = getPnFromHeight(pt, t, b, h, h_max, n_cells, swr, snr, pW, pCap, g, rhow, rhog, varargin)
    % Internal function: Reconstruct water pressure from VE assumption
    % T is top of column for each cell, t is top of cells to be
    % reconstructed, with b/B having the same interpretation for the
    % respective cell bottoms. h is the height of the gas plume in each
    % cell.
    % pW is UPSCALED water pressure at the bottom of the cell,
    % pCap is PSEUDO capillary pressure: pCap = pCap(h, H)
    % g the gravity magnitude and rhow/rhog the water and gas densities.

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
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0, 'transition', [], 'cNVE', []);
    opt = merge_options(opt, varargin{:});
       
    dz_above = max(pt - (t+h_max), 0);     
    h_above = min(pt - t, h);
    hmax_above = (pt - t) - h_above - dz_above;
              
    if strcmp(opt.transition, 'horizontal') || strcmp(opt.transition, 've') || strcmp(opt.transition, 'vertical') % saturation-weighted for horizontal ve-to-fine   
        % Assuming non-zero phase pressure in residual zone
        pG = pW + pCap + rhow.*g.*(swr.*h_above + (1-snr).*hmax_above + dz_above) ...
                        + rhog.*g.*((1-swr).*h_above + snr.*hmax_above);
                    
        % Odd's suggestion: no phase pressure in residual zone
%         pG = pW + pCap + rhow.*g.*((1-snr).*hmax_above + dz_above) ...
%                         + rhog.*g.*(1-swr).*h_above; 
        
        if any(opt.cNVE)
           % Special treatment for NVEs where residual sat filled from below
           cNVE = ismember(n_cells, opt.cNVE);

           pt_t = pt(cNVE) - t(cNVE); % height from top to point pt
           hmax_above = max(pt(cNVE) - (b(cNVE)-h_max(cNVE)), 0); % max(x, 0) to ensure no negatives
           dz_above = pt_t - hmax_above; 

           % Gas pressure in residual zone
           pG(cNVE) = pW(cNVE) + pCap(cNVE) + rhog(cNVE).*g.*snr.*hmax_above ...
                                            + rhow(cNVE).*g.*((1-snr).*hmax_above + dz_above);
           % No gas pressure in residual zone
           %pG(cNVE) = pW(cNVE) + pCap(cNVE) + rhow(cNVE).*g.*pt_t; % only water is mobile in cNVE cells --> no pressure contribution from hydrostatic co2 profile
        end
    else
        %pG = pW + pCap + rhog.*g.*(pt - t); % dupuit assumption for all other transitions (+ between pCap and rhog because upscaled gas pressure is defined at top)
        pG = pW + pCap; % transition between same discretization -> use original formula
    end
  
end
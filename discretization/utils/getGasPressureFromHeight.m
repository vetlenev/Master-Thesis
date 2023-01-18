function pG = getGasPressureFromHeight(pt, t, b, h, h_max, pW, pCap, g, rhow, rhog, varargin)
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
    opt = struct('ve_model', 'sharp_interface', 'p_entry', 0, 'transition', []);
    opt = merge_options(opt, varargin{:});
    
    %dz_below = (b - pt); % part of column between bottom of virtual cell and point of interest pt          
    %h_below = max(t + h - pt, 0); % part of plume below point pt in virtual cell   
    dz_above = max(pt - (t+h), 0);     
    h_above = min(pt - t, h);
              
    if strcmp(opt.transition, 'horizontal') || strcmp(opt.transition, 've') || strcmp(opt.transition, 'vertical') % saturation-weighted for horizontal ve-to-fine
        %pW_local = pW - ((dz_below - h_below).*rhow + h_below.*rhog).*g + opt.p_entry; 
        
        %pG = pW + pCap - ((dz_below - h_below).*rhow + h_below.*rhog).*g;
        pG = pW + pCap + (rhow.*dz_above + rhog.*h_above).*g;
%         if any(h_below > 0)
%             pG(h_below > 0) = pG(h_below > 0) - opt.p_entry;
%         end
        
    else
        %pG = pW + pCap + rhog.*g.*(pt - t); % dupuit assumption for all other transitions (+ between pCap and rhog because upscaled gas pressure is defined at top)
        pG = pW + pCap; % transition between same discretization -> use original formula
    end

end
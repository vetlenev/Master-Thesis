function [max_diff, max_diff_year, max_CO2] = plotRVEDifference(ax, vols_hybrid, vols_fine, t, discr_region, varargin)
% Generate a trapping inventory plot from a simulation result.  
% 
% The simulation result (set of states) first needs to be repackaged using the
% 'makeReports' function.
%
% SYNOPSIS:
%   function plotTrappingDistribution(ax, report, varargin)
%
% DESCRIPTION:
%
% PARAMETERS:
%   ax       - handle to figure into which to draw the plot
%   report   - simulation results, repackaged using the 'makeReports' function
%   varargin - optional argument pairs allowing to specify the location and
%              orientation of the legend. 
%
% SEE ALSO:
%   `makeReports`.

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

    opt.legend_location = 'northwest';
    opt.legend_orientation = 'vertical';
    opt.logScale = false;
    opt = merge_options(opt, varargin{:});
    
    % legend names
    names = discr_region;
    
    % Permuting volumes to comply with the ordering we want for reporting
    vols_diff = abs(vols_hybrid - vols_fine);       

    % Max difference between models through simulation
    [max_diff, max_idx] = max(sum(vols_diff, 2));
    max_CO2 = sum(vols_fine, 2);
    max_CO2 = max_CO2(max_idx);
    max_diff_year = round(t(max_idx));
    
    % Plotting trapping history
    area(ax, t, vols_diff./max_diff*100);  
    
    % Setting consistent colors
    col = getInventoryColors([7 6 5 5 4 3 2 1]);
    col(4,:) = 1/2 * (col(3,:) + col(5,:));
    chld = get(gca, 'Children');
    cur_col = 1;
%     for i = 1:numel(chld)
%         obj = get(chld(i));
%         % R2014 and later releases use Type=area, whereas earlier releases
%         % use Type=hggroup.
%         if strcmpi(obj.Type, 'hggroup') || strcmpi(obj.Type, 'area')
%             set(chld(i), 'FaceColor', col(cur_col, :));
%             
%             if (cur_col == 4) && skip_subscale
%                 set(get(get(chld(i), 'Annotation'), 'LegendInformation'), ...
%                     'IconDisplayStyle', 'off');
%             elseif (cur_col == 8) && skip_dissolution
%                 set(get(get(chld(i), 'Annotation'), 'LegendInformation'), ...
%                     'IconDisplayStyle', 'off');
%             end
%             cur_col = cur_col + 1;
%         end
%     end
        
    % Other formatting issues   
    axis tight
    if opt.logScale
        set(gca, 'XScale', 'log')
        xlim([0.1, max(t)])
    end
    xlabel('Years since simulation start');
    ylabel('Percentage')
    if numel(names) <= 10
        legend(names, ...
           'location', opt.legend_location, ...
           'orientation', opt.legend_orientation);
    end
end

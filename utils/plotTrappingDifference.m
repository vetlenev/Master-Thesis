function [max_diff, max_diff_year] = plotTrappingDifference(ax, report_hybrid, report_fine, varargin)
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

    opt.legend_location = 'eastoutside';
    opt.legend_orientation = 'vertical';
    opt.logScale = false;
    opt = merge_options(opt, varargin{:});
    
    % legend names
    names = {'Dissolved'           , 'Structural residual' , ...
             'Residual'            , 'Residual in plume'  , ...
             'Structural subscale' , 'Structural plume'   , ...
             'Free plume'          , 'Exited'};
    
        
    % Loading and extracting trapping history
    t_hist = [report_hybrid.t]';
    t_hist_ = [report_fine.t]';
    assert(isequal(t_hist, t_hist_));

    mass_fine = reshape([report_fine.masses]',8, [])';  
    mass_hybrid = reshape([report_hybrid.masses]',8, [])'; 
    mass_diff = convertTo(abs(mass_fine - mass_hybrid), mega*kilo);

    [max_diff, max_diff_year] = max(sum(mass_diff, 2));
    max_diff_year = t_hist(max_diff_year)/year();
    mass_diff = mass_diff./max_diff*100;
    
    % Permuting volumes to comply with the ordering we want for reporting
    P = [1 2 3 4 6 5 7 8];
    mass_diff = mass_diff(:,P);

    % Checking if subscale trapping and dissolution is present
    sum_tmp = sum(mass_diff);
    skip_subscale = (sum_tmp(5) == 0);
    skip_dissolution = (sum_tmp(1) == 0);
    if skip_subscale
        names(5) = [];
    end
    if skip_dissolution
        names(1) = [];
    end
    
    % Plotting trapping history
    area(ax, ...
         convertTo(t_hist, year), ...
         [mass_diff(:, 1:7), mass_diff(:,8)]);
         % convertTo([mass_hist(:, 1:7), mass_hist(:,8)-sum(mass_hist(:,1:7),2)], ...
         %           mega*kilo));
    
    % Setting consistent colors
    col = getInventoryColors([7 6 5 5 4 3 2 1]);
    col(4,:) = 1/2 * (col(3,:) + col(5,:));
    chld = get(gca, 'Children');
    cur_col = 1;
    for i = 1:numel(chld)
        obj = get(chld(i));
        % R2014 and later releases use Type=area, whereas earlier releases
        % use Type=hggroup.
        if strcmpi(obj.Type, 'hggroup') || strcmpi(obj.Type, 'area')
            set(chld(i), 'FaceColor', col(cur_col, :));
            
            if (cur_col == 4) && skip_subscale
                set(get(get(chld(i), 'Annotation'), 'LegendInformation'), ...
                    'IconDisplayStyle', 'off');
            elseif (cur_col == 8) && skip_dissolution
                set(get(get(chld(i), 'Annotation'), 'LegendInformation'), ...
                    'IconDisplayStyle', 'off');
            end
            cur_col = cur_col + 1;
        end
    end
        
    % Other formatting issues   
    axis tight
    if opt.logScale
        set(gca, 'XScale', 'log')
        xlim([0.1, max(convertTo(t_hist, year))])
    end
    xlabel('Years since simulation start');
    ylabel('Percentage (%)')
    legend(names, ...
           'location', opt.legend_location, ...
           'orientation', opt.legend_orientation);
end

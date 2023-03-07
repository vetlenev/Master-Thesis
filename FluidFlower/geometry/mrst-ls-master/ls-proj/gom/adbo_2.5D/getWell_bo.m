function [W, timesteps, wellInx] = getWell_bo(G, rock, fluid, wells, t)
%
%
%

% Individualize vars
injTime = t(1);      simTime = t(2);
mrate   = wells.mrate; nwell = wells.num;

% Well(s)
wells.loc(1) = max(G.faces.centroids(:, 1))/2;
wellLoc = wells.loc;
dirw = 'z';
dist = pdist2(G.cells.centroids, wellLoc);
[~, wellInx] = min(dist,[],1);

injmass = mrate*injTime;
rhoInj = fluid.rhoGS;
injrate = injmass/(rhoInj*t(1)*nwell);                                      % [Sm^3/s]
W = addWell([ ], G, rock, wellInx, 'Name', 'I1', 'Dir', dirw, ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1], ...            % order is always 'WOG' ('OG' in GenericBlackOilModel)
            'refDepth', G.cells.centroids(wellInx, G.griddim));             % refDepth must be included for 'bhp' wells.

% Timesteps
% if isfield(fluid, 'krHyst') && fluid.krHyst ~= 0 || ...
%         isfield(fluid, 'pcHyst') && fluid.pcHyst ~= 0
%     if injTime/year == 20 && simTime/year == 100
%         %reportTimes = [1*hour, [1, 2, 7, 14, 21, 30, 60, 90, 120, 150, 180, ...
%         %    240, 300, 365]*day, [1+(1/12):1/12:10, 10+(2/12):2/12:23, ...
%         %    23.25:0.25:30 30.5:0.5:60, 61:100]*year];
%         reportTimes = [[1, 2, 7, 14, 21, 30, 60, 90, 120, 150, 180, ...
%             240, 300, 365, 456.25, 547.5, 638.75, 730, 821.25]*day, ...
%             [2.5:0.5:22.5, 23:40, 42:2:60, 65:5:100]*year];
%     end
%     
% else
    if injTime/year == 30 && simTime/year == 200
        reportTimes = [1*hour, [1, 2, 7, 14, 21, 30, 60, 90, 120, 150, 180, ...
             240, 300, 365, 456.25, 547.5, 638.75, 730, 821.25]*day, ...
            [2.5:0.5:10, 11:28, 28.5:0.5:32, 33:1:50, 52:2:70 75:5:150 160:10:200]*year];
    elseif injTime/year == 50 && simTime/year == 500
        reportTimes = [1*hour, [1, 2, 7, 14, 21, 30, 60, 90, 120, 150, 180, ...
             240, 300, 365, 456.25, 547.5, 638.75, 730, 821.25]*day, ...
            [2.5:0.5:10, 11:48, 48.5:0.5:52, 53:1:60, 62:2:100 105:5:200 210:10:500]*year];
    end
%end
timesteps = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==simTime, 'sum of timesteps must equal simTime')

end
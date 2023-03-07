function [sPoro, cPoro] = endMemberPoro(G, rock, ucids, varargin)
%
%
%


% Get variable inputs
opt = struct('sandOpt', 'constant', ...
             'clayMineral', 'kaolinite', ...
             'rhoBrine', 1050); 
opt = merge_options(opt,varargin{:});

% Pure sand (Note that this neglects cataclasis in the sand).
% Mean porosity of each sand unit is assigned to mean depth of each sand 
% unit. Then pure-sand porosity in the fault is linearly interpolated from 
% these values, based on depth.
spor = 0.4;
ssid = ucids.Vcl <= spor;
disp(['endMemberPoro.m: Linear interpolation for fault poro considers Vcl <= ' num2str(spor) ' to be sand.'])
depth = zeros(numel(ssid), 1);
poroAvg = depth;
for n = 1:numel(ucids.Vcl)
        depth(n) = mean(G.cells.centroids(ucids.unit_cell_ids{n}, G.griddim));
        poroAvg(n) = mean(rock.poro(ucids.unit_cell_ids{n}));
end
poroAvg(~ssid) = nan;

switch opt.sandOpt
    case 'constant'
        if sum(~isnan(poroAvg)) == 1 % cannot interpolate
            depth = [2.4982; 1.0360]*10^3;
            ssid  = 1:2;
            poroAvg = [0.2655; 0.2799];
        end
        poroAvg = interp1(depth(ssid), poroAvg(ssid), depth, 'linear');
        if isfield(ucids, 'upper') && ucids.upper == 2
            poroAvgf = repelem(poroAvg, 2);
            depthf   = [max(G.cells.centroids([ucids.unit_cell_ids{3:4}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{3:4}], 3));
                        max(G.cells.centroids([ucids.unit_cell_ids{5}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{5}], 3))];
        elseif isfield(ucids, 'upper') && ucids.upper == 1
            poroAvgf = repelem(poroAvg, 2);
            depthf = [max(G.cells.centroids([ucids.unit_cell_ids{5:6}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{5:6}], 3));
                      max(G.cells.centroids([ucids.unit_cell_ids{7}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{7}], 3));
                      max(G.cells.centroids([ucids.unit_cell_ids{8:9}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{8:9}], 3));
                      max(G.cells.centroids([ucids.unit_cell_ids{10}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{10}], 3))];
        else
            poroAvgf = repelem(poroAvg(4:end), 2);
            depthf = [max(G.cells.centroids([ucids.unit_cell_ids{8:9}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{8:9}], 3));
                      max(G.cells.centroids([ucids.unit_cell_ids{10}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{10}], 3));
                      max(G.cells.centroids([ucids.unit_cell_ids{11:12}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{11:12}], 3));
                      max(G.cells.centroids([ucids.unit_cell_ids{13}], 3)); min(G.cells.centroids([ucids.unit_cell_ids{13}], 3))];
        end
        [xData, yData] = prepareCurveData( depthf, poroAvgf );
        sPoro = fit(xData, yData, 'linearinterp');                          % sand porosity (fcn of depth)
        
    case 'linearInterp'
        [xData, yData] = prepareCurveData( depth(ssid), poroAvg(ssid) );
        sPoro = fit(xData, yData, 'linearinterp');                          % sand porosity (fcn of depth)
end

% Pure clay (based on compaction curves since clay poro much more dependent
% on compaction than shear strain). Compaction data and equations in 
% Revil et al. (2002) are used.
rhog = 2650;                % bulk density of the grains
rhob = opt.rhoBrine;        % bulk density of the fluid
switch opt.clayMineral
    case 'kaolinite'
        phi0 = 0.55;                                                        % depositional porosity
        bm   = 6*10^(-8);                                                   % 1/Pa "compaction coefficient of the shale end-member"
        zm   = 1/(norm(gravity)*bm*(rhog-rhob));                            % characteristic length [m]
end
cPoro = @(z) phi0./(phi0 + (1-phi0)*exp(z./zm));                            % clay porosity (fcn of depth)

end
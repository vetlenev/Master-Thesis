function [W, timesteps, wellInx] = getWells(G, G_dat, rock, opt, model_case, inj_type)
%
%
%

% wellInx
idwG = find(~isnan(G_dat.wellNo));
ninj = size(opt.inj_loc, 1);
wellInx = zeros(1, ninj);
W = [];
for n=1:size(opt.inj_loc, 1)
    dist = pdist2(G.cells.centroids(idwG,2:3), opt.inj_loc(n,:));
    [~, inx] = min(dist,[],1);
    wellInx(n) = idwG(inx);
    
    % Define wells
    injrate = opt.rate{(n-1)*2+1,1}(end)*(milli*litre)/(minute); 
    W = addWell(W, G, rock, wellInx(n), 'Name', ['I' num2str(n)], ...
                'Dir', 'z', 'Type', 'rate', 'Val', injrate, ...
                'compi', [0, 1], ...        % order is 'OG'
                'refDepth', G.cells.centroids(wellInx(n), G.griddim), ...
                'Radius', 0.9*1e-3);              
end

% Other quantities (later)
%injvol = injrate*t(2);
%rhoInj = fluid.rhoGS;
%mrate = injvol*rhoInj/(t(2)*wellno);   

% Report times
tsim = opt.schedule{end}*minute;
if model_case < 4
    if inj_type == 1
        % Get variables
        rd1.end = opt.schedule{2};
        rd1.t = max(opt.rate{2,2});
        rd1.tstep = opt.rate{2,2};
        i2_t0 = opt.schedule{3}(1)*minute + opt.schedule{3}(2);
        i2_t0min = i2_t0/minute;
        i2_tru_end = i2_t0min+max(opt.rate{3,2});   % ramp up end
        rd2.end = opt.schedule{4};
        rd2.t = sum(opt.rate{4,2})/minute;
        rd2.tstep = opt.rate{4,2};
        rd2.start = rd2.end - rd2.t;
        rd2.startInt = fix(rd2.start);
        rd2.startSep = [rd2.startInt (rd2.start-rd2.startInt)*minute];
        rd2.end = ((rd2.start*minute)+sum(rd2.tstep))/minute;
        
        trep = [(1:10)*minute, ...                             % rampup inj rate + 4min
            (11:1:rd1.end-rd1.t)*minute, ...               % inj 1
            ((rd1.end-rd1.t)+rd1.tstep)*minute, ...        % rampdown inj 1
            [(opt.schedule{2}+1:opt.schedule{3}(1)-1)*minute i2_t0], ...    % between I1 and I2
            (i2_t0min+1:i2_tru_end+5)*minute, ...          % rampup inj 2
            (i2_tru_end+6:1:rd2.start)*minute, ...         % inj 2
            (rd2.start*minute)+cumsum(rd2.tstep), ...      % rampdown inj 2
            [rd2.end+6:5:550 560:10:opt.schedule{end}]*minute];

    elseif inj_type == 2
        rd1.end = opt.schedule{2};
        rd1.t = sum(opt.rate{2,2});
        rd1.d1 = [283 44];
        rd1.tstep = opt.rate{2,2};
        
    trep = [(1:10)*minute, ...                             % rampup inj rate + 4min
            [(11:1:rd1.d1(1))*minute rd1.d1(1)*minute + rd1.d1(2)], ...     % inj 1
            rd1.d1(1)*minute + rd1.d1(2) + cumsum(rd1.tstep), ... % rampdown
            [rd1.end(1)+1:1:300 302:2:350 355:5:600 ... % rest
            610:5:opt.schedule{end}]*minute];
    elseif inj_type == 3   
        % Get variables
        rd1.end = opt.schedule{2};
        rd1.t = max(opt.rate{2,2});
        rd1.tstep = opt.rate{2,2};
        rd1.start(1) = rd1.end(1) -rd1.t;
        rd1.start(2) = rd1.end(2);
        i2_t0 = opt.schedule{3}(1)*minute + opt.schedule{3}(2);
        i2_t0min = i2_t0/minute;
        i2_tru_end = i2_t0min+max(opt.rate{3,2});   % ramp up end
        rd2.end = opt.schedule{4};
        rd2.t = max(opt.rate{2,2});
        rd2.tstep = opt.rate{4,2};
        rd2.start(1) = rd2.end(1) - rd2.t;
        rd2.start(2) = rd2.end(2);
        
        trep = [(1:10)*minute, ...                                                   % rampup inj rate + 4min
                [(11:1:rd1.start(1))*minute rd1.start(1)*minute + rd1.start(2)], ... % inj 1
                (rd1.start(1)+rd1.tstep)*minute + rd1.start(2), ...                  % rampdown inj 1
                [(rd1.end(1)+1:opt.schedule{3}(1))*minute i2_t0], ...                % between I1 and I2
                 (i2_t0min+1:i2_tru_end+5)*minute, ...                               % rampup inj 2
                [(i2_tru_end+6:1:rd2.start(1))*minute rd2.start(1)*minute + rd2.start(2)], ...  % inj 2
               (rd2.start(1)+rd2.tstep)*minute + rd2.start(2), ...                   % rampdown inj 2
                [rd2.end(1)+4:5:550 560:10:opt.schedule{end}]*minute];
    end
                           
timesteps = [trep(1) diff(trep)];
assert(sum(timesteps)== tsim, 'sum of timesteps must equal simTime')
end

end
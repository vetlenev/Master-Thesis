function [W, timesteps, wellInx] = bWells(inj_case, G, G_dat, rock, opt)
%
%
%

% wellInx
idwG = find(~isnan(G_dat.wellNo));
nports = opt.wellno;
if inj_case < 3  % well test 1 or 2
    if inj_case == 1
        wellInx = idwG(4);          % (17,7), (9, 3) is 1.
    elseif inj_case == 2
        wellInx = idwG(1);   
    end
    injrate = opt.rate{1,1}(end)*(milli*litre)/(minute); 
    W = addWell([], G, rock, wellInx, 'Name', 'I1', ...
                        'Dir', 'x', 'Type', 'rate', 'Val', injrate, ...
                        'compi', [1, 0], ...        % order is 'OG'
                        'refDepth', G.cells.centroids(wellInx, G.griddim), ...
                        'Radius', 0.9*1e-3); 
            
elseif inj_case == 3    % CO2 inj
    wellInx = zeros(1, nports);
    W = [];
    for n=1:nports
        dist = pdist2(G.cells.centroids(idwG,2:3), opt.injloc(n,2:3));
        [~, inx] = min(dist,[],1);
        wellInx(n) = idwG(inx);
        injrate = opt.rate{(n-1)*2+1,1}(end);
        %injrate = opt.rate{(n-1)*2+1,1}(end)*(milli*litre)/(minute);
        W = addWell(W, G, rock, wellInx(n), 'Name', ['I' num2str(n)], ...
                    'Dir', 'x', 'Type', 'rate', 'Val', injrate, ...
                    'compi', [0, 1], ...        % order is 'OG'
                    'refDepth', G.cells.centroids(wellInx(n), G.griddim), ...
                    'Radius', 0.9*1e-3);
    end
end

% Other quantities (later)
%injvol = injrate*t(2);
%rhoInj = fluid.rhoGS;
%mrate = injvol*rhoInj/(t(2)*wellno);   

% Report times
tsim = opt.schedule(end)*minute;
if inj_case == 1 || inj_case == 2
    ti0 = opt.schedule([1 3 5]);
    tif = opt.schedule([2 4 6]);
    trep = [opt.rate{1,2} ([1 2 5:5:tif(1)])*minute ...
            tif(1)*minute+opt.rate{2,2} ([tif(1)+1 tif(1)+5:5:ti0(2)])*minute ...
            ti0(2)*minute+opt.rate{3,2} ([ti0(2)+1 ti0(2)+5:5:tif(2)])*minute ...
            tif(2)*minute+opt.rate{4,2} ([tif(2)+1 tif(2)+5:5:ti0(3)])*minute ...
            ti0(3)*minute+opt.rate{5,2} ([ti0(3)+1 ti0(3)+5:5:tif(3)])*minute ...
            tif(3)*minute+opt.rate{6,2} ([tif(3)+1 tif(3)+5:5:opt.schedule(end)])*minute];
    
elseif inj_case == 3
    if strcmp(opt.rampdown, 'fast')
       trep = [opt.rate{1, 2} ...                        % rampup inj 1
               (1:1:135)*minute ...                      % inj 1 only
               135*minute + opt.rate{3,2} ...            % inj 2 rampup
               (136:1:300)*minute ...                    % inj 1 and inj 2
               300*minute + opt.rate{4,2} ...            % inj 1 and inj 2 rampdown
               [301 (302:2:360)]*minute ...              % 2 min until 6h
               (365:5:720)*minute ...                    % 5 min until 12h
               (730:10:5*24*60)*minute];                 % 10 min unitl 5 days
    else
       trep = [opt.rate{1, 2} ...                        % rampup inj 1
               (1:1:135)*minute ...                      % inj 1 only
               135*minute + opt.rate{3,2} ...            % inj 2 rampup
               (136:1:298)*minute ...                    % inj 1 and inj 2
               (298 + opt.rate{4,2})*minute ...          % inj 1 and inj 2 rampdown
               (314:2:360)*minute ...                    % 2 min until 6h
               (365:5:720)*minute ...                    % 5 min until 12h
               (730:10:5*24*60)*minute];                 % 10 min unitl 5 days
    end
end

timesteps = [trep(1) diff(trep)];
assert(sum(timesteps)== tsim, 'sum of timesteps must equal simTime')
end
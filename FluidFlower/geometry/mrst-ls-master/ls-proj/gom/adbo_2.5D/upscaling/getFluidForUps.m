function [rock, fluid, fault] = getFluidForUps(G, rock, fault, opt)
%
%
%
pth = fullfile(mrstPath('ls-proj'), 'gom/adbo_2.5D/eclipse_data_files');
rock.regions.rocknum = ones(G.cells.num, 1);

if strcmp(opt.fault, 'test')
    rock.regions.saturation = ones(G.cells.num, 1);  % sand
    permvals = unique(rock.perm);
    id_lowperm = rock.perm == min(permvals);
    rock.regions.saturation(id_lowperm) = 2;         % clay
    rs = rock.regions.saturation;
    fault.Grid.isSmear = false(G.cells.num, 1);
    fault.Grid.isSmear(rs == 2) = true;
    
    % load deck and initialize fluid
    %fn  = fullfile(pth, 'gom_co2brine_2reg_synt.DATA');
    fn  = fullfile(pth, 'gom_forUps_theta30.DATA');
    deck = convertDeckUnits(readEclipseDeck(fn));
    deck.REGIONS.ROCKNUM = rock.regions.rocknum;
    fluid = initDeckADIFluid(deck);
    fluid.isclay = [0 1];
    
elseif strcmp(opt.strati, 'GoM')
    %domain = flipud(fault.MatMap.units)';
    %domain = domain(:);
    if strcmp(opt.fault, 'predict')
        isclay = flipud(fault.MatMap.vals)'; 
        isclay = isclay(:);
    elseif strcmp(opt.fault, 'predict_3D')
        isclay = double(fault.Grid.isSmear);
    end

    % Find rock.regions.saturation, accounting for clay or sand in each
    % domain (number of regs >= MatMap.units). Note that PREDICT only
    % outputs in fault.MatMap.units the sand (sand domain) or clay (clay
    % domain).
%    rock.regions.saturation = zeros(G.cells.num, 1);
%     for j = 1:numel(fault.MatMap.unit)
%         idu = domain == fault.MatMap.unit(j);
%         idc = all([idu, isclay], 2);
%         unitInGapsColum = repelem(fault.MatMap.unitInClayGaps(j), ...
%                                   G.cartDims(1), 1);
%         unitInClayGaps = zeros(G.cartDims(1));
%         unitInClayGaps = full(spdiags(repmat(unitInGapsColum, 1, fault.MatMap.nDiag(j)), ...
%                                       fault.MatMap.DiagBot(j):fault.MatMap.DiagTop(j), ...
%                                       unitInClayGaps));
%         unitInClayGaps = flipud(unitInClayGaps)';
%         unitInClayGaps = unitInClayGaps(:); % same-clay with different sand in gaps
%         ids = all([idu, ~isclay, unitInClayGaps], 2);
%         unitEntries = sum(fault.MatMap.unit == fault.MatMap.unit(j));
%         if sum(idc) > 0 && sum(ids) > 0  || ...
%            unitEntries > 1 && fault.MatMap.isclay(j)     % smear and sand
%             rock.regions.saturation(idc) = fault.MatMap.unit(j);
%             rock.regions.saturation(ids) = fault.MatMap.unitInClayGaps(j);
%         elseif sum(idc) > 0  || sum(ids) > 0    % only smear or only sand
%             rock.regions.saturation(idu) = fault.MatMap.unit(j);
%         end
%     end
    rock.regions.saturation = fault.Grid.units;
    rock.regions.rocknum(fault.Grid.isSmear) = 2;
    
    fn  = fullfile(pth, 'gom_forUps_theta30_PVDO_incompRock.DATA');
    deck = convertDeckUnits(readEclipseDeck(fn));
    deck.REGIONS.ROCKNUM = rock.regions.rocknum;
    fluid = initDeckADIFluid(deck);
    
    % Get pore compressibility
    cP = deck.PROPS.ROCK(:, 2);
    fluid = assignPvMult(fluid, cP, rock.regions.rocknum);
    
    if opt.theta(2) == 60       % after kr for theta 60 finished, substitute for .DATA file
        f = fluid;  fluid = [];
        sg_clay = linspace(0, 0.75, 30);
        sg_clay = [sg_clay(1) 1e-3 sg_clay(2:end-1) 0.75-0.001 sg_clay(end)];
        pc_smooth = @(x) 6.265e+6*x.^2.153 + 162379.542066555;
        pc_clay = [0 pc_smooth(sg_clay(2:end))];
        fluid.pcOG = {f.pcOG{1},  @(sg) interp1(sg_clay, pc_clay, sg)};
        fluid.krPts = f.krPts;
        fluid.krPts.g(2,3) = 1-0.25; fluid.krPts.og(2,2) = 0.25;
    end
    
    %Scale Pc for each unit (we neglect the small variation introduced in
    %PREDICT in poro and perm for cells in same unit)
    id_reg = unique(rock.regions.saturation);
    if strcmp(opt.dir, 'x'), id_perm = 1;    % xx
    else, id_perm = numel(rock.perm(1,:));                       % zz
    end
    refPerm = [7.60393535652603e-13, 0.00208*milli*darcy];   % [sand, clay]
    refPoro = [0.289875, 0.07];         % [    "     ]
    pce_rmse = 0.2953;
    % pce_hg: sample randomly (uniform law) between +- rmse
    pce_hg = @(log10_k_md, r) 10.^(-0.1992*log10_k_md + 1.407 ...
                                    -pce_rmse + r*2*pce_rmse);  % bar, pce_model_clay.m
    pce_co2wat = @(pce_hg, theta) 1e5*pce_hg*abs(cosd(theta)*25/(cosd(140)*485)); % Pa
    nreg = numel(id_reg);
    pcOGups = cell(1, max(id_reg));
    fluid.isclay = false(1, nreg);
    for n=1:nreg
        id_cells = rock.regions.saturation == id_reg(n);
        permv = mean(rock.perm(id_cells, id_perm));     % m2
        porov = mean(rock.poro(id_cells));
        if unique(isclay(id_cells))         % is clay smear
            log10_k_md = log10(permv/(milli*darcy));
            refPc = fluid.pcOG{2};
            fluid.isclay(n) = true;
            rv = rand(1);
            pcOGups{id_reg(n)} = @(sg) refPc(sg) * ...
                  pce_co2wat(pce_hg(log10_k_md, rv), opt.theta(2))/refPc(1e-3);
        elseif ~unique(isclay(id_cells))    % is sand
            refPermv = refPerm(1);
            refPorov = refPoro(1);
            refPc = fluid.pcOG{1};
            pcOGups{id_reg(n)} = @(sg) refPc(sg).*sqrt((refPermv*porov)./(refPorov*permv));
        end
        
    end
    fluid.pcOG = pcOGups;
end

end
function [pc, sg, fluid, sgmax, pc_fcn, pc_max, pc_min, sat_min, sat_max] = ...
                    upscalePcReg(G, fluid, rock, opt, makeplots)
%
% For CL method: based on upscalePCCaplimit, includes support for multiple
% and a fluid initialized with initDeckADIFluid. Assumptions:
%  1. The two phases are O and G (so that fluid.pcOG is used).
%  2. the fluid was obtained using initDeckADIFluid and a .DATA input file
%     with SGOF kr-Pc data.
%
% If need to add gravity, follow upPcOW.m in the steady-state module. For
% now, no gravity support (the upscaling process neglects gravity).
%
% Inputs:
%   - G: MRST grid structure
%   - fluid: MRST fluid structure (from .DATA input file with SGOF kr-Pc
%                                  data)
%   - reg: rock.regions.saturation field (Nc x 1). Note this is for
%          drainage only (assumes reversible Pc curves).
%
% Outputs:
%   - See code
%
   
    % Number of regions and cells in region
    reg = rock.regions.saturation;
    id_reg = unique(rock.regions.saturation);
    nreg = numel(id_reg);
    ncreg = sum(reg(:)==(1:max(reg)));
    ncreg(ncreg==0) = [];
    
    %
    fluid.pcInv = cell(1, max(reg));
    pc_min = 10^10;
    pc_max = 0;
    %pce_max = 0;
    pcv2 = 1e+6;
    sgmax = zeros(1, nreg);
    for n=1:nreg
        % Add inverse pc field
        if isempty(opt.sg)
            sgmin = fluid.krPts.g(n,1);
            sgmax(n) = fluid.krPts.g(n,3);
            sgvals = linspace(sgmin, sgmax(n), pow2(6)-1)';
        elseif strcmp(opt.sg, 'sandClay')
            if fluid.isclay(n) == 1
                sgmin = fluid.krPts.g(2,1);
                sgmax(n) = fluid.krPts.g(2,3);
                sgvals = linspace(sgmin, sgmax(n), pow2(6)-1)';
            else
                sgmin = fluid.krPts.g(1,1);
                sgmax(n) = fluid.krPts.g(1,3);
                sgvals = linspace(sgmin, sgmax(n), pow2(6)-1)';
            end
            
        end
        sgvals = [sgvals(1); sgvals(1)+1e-3; sgvals(2:end)];
        pcvals = fluid.pcOG{id_reg(n)}(sgvals);
        pcv2 = min([pcv2, pcvals(2)]);
        fluid.pcInv{id_reg(n)} = @(pcOG) interp1(pcvals, sgvals, pcOG);
        
        % hack to find maximum/minimum capillary pressure
        % assumes finescale pc gives the maximum/minimum out of range
        pc_min = min([pc_min fluid.pcOG{id_reg(n)}(sgmin)]);
        pc_max = max([pc_max fluid.pcOG{id_reg(n)}(sgmax(n))]);
        %pce_max = max([pce_max fluid.pcOG{n}(1e-3)]); %  assumes pce is at sg=1e-3 in .DATA file
    end
    
    if strcmp('oper', opt.pc_mode)  
        % ----------------- ordinary percolation --------------------------
        %pc_val = linspace(pc_min, pc_max, pow2(8)-1);
        %diffv = pc_val(2)-pc_val(1);
        %pc_val = [pc_val(1) min([pc_val(1)+0.1*diffv pcv2]) pc_val(2:end)];
        pc_val = logspace(log10(pcv2), log10(0.99*pc_max), pow2(8)-1);
        pc_val = [0 0.98*pcv2 pc_val];
    
        volume=sum(G.cells.volumes.*rock.poro);
        s_val=nan(numel(pc_val), 1);
        for i=1:numel(pc_val)
            s_cells = nan(G.cells.num, 1);
            for n=1:nreg
                idcell = reg == id_reg(n);
                if pc_val(i) > fluid.pcOG{id_reg(n)}(sgmax(n))     % pcInv returns NaN
                   pcv = fluid.pcOG{id_reg(n)}(sgmax(n));
                else
                   pcv = pc_val(i);
                end
                s_cells(idcell)=fluid.pcInv{id_reg(n)}(pcv*ones(ncreg(n), 1));
            end
            %idv = ~isnan(s_cells);
            %volume=sum(G.cells.volumes(idv).*rock.poro(idv));
            s_val(i)=sum(s_cells.*G.cells.volumes.*rock.poro, 'omitnan')/volume;
        end
        s_pce = 1e-5;
        
    elseif strcmp('inv-per', opt.pc_mode)
        %
        % ---------------- Invasion-percolation ---------------------------
        %
        pc_val = logspace(log10(pcv2), log10(0.99*pc_max), pow2(6)-2);
        pc_val = [0 0.98*pcv2 pc_val];
        t = opt.t;
        
        % Entry pressures for all cells
        pce = zeros(G.cells.num, 1);
        for n=1:nreg
            pce(rock.regions.saturation == id_reg(n)) = fluid.pcOG{id_reg(n)}(1e-3);
        end
        
        % Set outer loop
        tic
        volume=sum(G.cells.volumes.*rock.poro);
        f2cn = gridCellNo(G);
        regs = unique(rock.regions.saturation);
        s_val = zeros(numel(pc_val), 1);
        sg = zeros(G.cells.num, 1);
        istrapped = false(G.cells.num, 1);
        idrem = false(numel(pc_val), 1);
        s_last = 0;
        maxIts = fix(5*G.cartDims(1));
        id_percolation = [];
        openCellsAll = [];
        for k=2:numel(pc_val)
        
            % Set inner loop
            evaluatedCells = [];
            pcv = pc_val(k);
            it = 1;
            while it < maxIts
                if it == 1
                    connectedCells = (1:G.cartDims(1)*G.cartDims(2))';
                else
                    %Here
                    connectedCells = G.faces.neighbors(G.cells.faces(ismember(f2cn,openCells)), :);
                    connectedCells = unique(reshape(connectedCells, numel(connectedCells), 1));
                    connectedCells = connectedCells(~ismember(connectedCells,[0; evaluatedCells]));
                end
                dpmax = max(pcv-pce(connectedCells));
                if dpmax <= 0
                    openCells = connectedCells(pcv >= pce(connectedCells));
                else
                    openCells = connectedCells(pcv >= pce(connectedCells) + (1-t)*dpmax);
                end 
                if isempty(openCells)
                    break
                elseif any(ismember(((G.cells.num-G.cartDims(1)*G.cartDims(2))+1):G.cells.num, ...
                                    openCells)) && isempty(id_percolation)
                       id_percolation = k;    
                end
                openCellsAll = [openCells; openCellsAll];
                for n=1:numel(openCells)    % needed bc trapping could vary around each cell
                    aroundCells = G.faces.neighbors(G.cells.faces(f2cn == openCells(n)), :);
                    aroundCells = unique(reshape(aroundCells, numel(aroundCells), 1));
                    aroundCells = aroundCells(~ismember(aroundCells, [0; openCells(n)]));
                    regid = rock.regions.saturation(openCells(n));
                    if sg(openCells(n)) < sgmax(regs==regid)-1e-3 && ~all(istrapped(aroundCells))
                        pcv_id = min(pcv, fluid.pcOG{regid}(sgmax(regs==regid)));
                        sg(openCells(n)) = fluid.pcInv{regid}(pcv_id);
                        if sg(openCells(n)) >= sgmax(regs==regid)-1e-3
                            istrapped(n) = true;
                        end
                    elseif ~istrapped(n)    % either surrounded or sgmax in cell
                        istrapped(n) = true;
                    end
                end
                evaluatedCells = unique([connectedCells; evaluatedCells]);
                it = it + 1;  
                if it == maxIts
                    error('maxIts reached but sweep not yet completed!')
                end
            end
            if isempty(id_percolation)
                s_val(k) = 0;           % sg=0 until percolation occurs.
            else
                s_val(k)=sum(sg.*G.cells.volumes.*rock.poro)/volume;
            end
            if abs(s_val(k)-s_last) < 1e-3      % only different s values
                idrem(k) = true;
            else
                s_last = s_val(k);
            end
%             figure(99); clf; plotGrid(G,'facecolor','none'); 
%             plotGrid(G, evaluatedCells, 'facecolor', 'r','edgecolor','none');
%             ax = gca;
%             ax.DataAspectRatio = [0.05 1 1];
%             ax.ZDir = 'normal';
%             view([30 20])

           if mod(k,5) == 0
              disp([num2str(k) '/' num2str(numel(pc_val)) ' completed.'])
           end
        end
        idrem(id_percolation-1) = false;
        idrem(id_percolation) = false;
        idrem(end) = false;
        s_pce = 1e-5;
        assert(s_pce < s_val(id_percolation));
        s_val(id_percolation-1) = s_pce;
        %pc_val(id_percolation-1) = pc_val(id_percolation);
        s_val(idrem) = [];
        pc_val(idrem) = [];
        toc
        
    elseif strcmp('sim_bo', opt.pc_mode) 
        % --------------- 2-phase, blackoil simulation (no g) -------------
        % This is way too slow to generate multiple curves, since it takes
        % several minutes to generate 1 point, even without gravity and
        % incompressible, immiscible fluids. 
        % TBD: Outer loop for all pc points, for now just tested with one.
        %
        pc_val = logspace(log10(pcv2), log10(0.99*pc_max), pow2(5)-1);
        pc_val = [0 pc_val];
        
        % We need to have same nreg for kr as for pc
        krg = fluid.krG;
        fluid.krG = cell(1, numel(fluid.isclay));
        fluid.krG(fluid.isclay) = krg(2);
        fluid.krG(~fluid.isclay) = krg(1);
        
        % We need kr, pc from 1 to end (otherwise mrst sim won't work).
        pcog = fluid.pcOG; 
        idu = find(~cellfun(@isempty, pcog));
        fluid.pcOG = cell(1, numel(idu));    fluid.pcOGu = pcog; 
        fluid.pcOG(1:numel(idu)) = pcog(idu);
        
        % And adjust rock.regions.saturation
        idv = 1:numel(idu);
        if sum(idv-idu) ~= 0
            for n=1:numel(idu)
                if idu(n) ~= idv(n)
                    rock.regions.saturation(rock.regions.saturation == idu(n)) = idv(n);
                end
            end
        end
        
        % Gravity
        gravity reset off
        
        
        % Initialization [introduce option for MRST grids to be generated
        % at zmax depth, and do porper initialization, for now we consider
        % it constant p.
        p_r = 9.8066*fluid.rhoOS*mean(opt.zmax{1}); % p at shallowest z
        p0 = repelem(p_r, G.cells.num, 1);
        s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
        rs0 = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
        rv0 = zeros(G.cells.num, 1);             % no vaporized water at the beginning
        state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
        
        
        % Model
        model = GenericBlackOilModel(G, rock, fluid, 'disgas', false, ...
                                    'vapoil', false, 'water', false);
                         
        % Acceleration and solver parameters
        model.AutoDiffBackend = DiagonalAutoDiffBackend('useMex', true, ...
                                              'deferredAssembly', true);
        model.toleranceCNV = 5e-3;
        model = model.validateModel();

        nls = getNonLinearSolver(model, 'TimestepStrategy', 'none', ...
                                 'useCPR', true);
        stepSel = StateChangeTimeStepSelector('targetProps', {'s'}, ...
                                              'targetChangeAbs', 0.2);
        nls.timeStepSelector = stepSel;
        nls.LinearSolver.maxIterations = 50;
        nls.useRelaxation = 1;
        nls.maxIterations = 8;
        nls.maxTimestepCuts = 12;
        nls.acceptanceFactor = 5;

        % Model changes (after acceleration!)
        model.operators.p0 = state0.pressure;                            % Needed for MyPvMult
        model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
        model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('RelativePermeability', MyRelativePermeability(model));
        model.FlowPropertyFunctions = model.FlowPropertyFunctions.setStateFunction('CapillaryPressure', MyBlackOilCapillaryPressure(model));
        model.minimumPressure = min(state0.pressure);

        
        % BCs
        f = false(G.faces.num, 1);
        f(boundaryFaces(G)) = true;
        f1 = all([G.faces.centroids(:, 2) < 1e-4, f], 2);
        topv = max(G.faces.centroids)-1e-4;
        f2 = all([G.faces.centroids(:, 2) > topv, f], 2);
        bc = addBC([], find(f1), 'pressure', p_r + mean(pc_val), 'sat', [0 1]);
        bc = addBC(bc, find(f2), 'pressure', 0, 'sat', [1 0]);
        %plotGrid(G); hold on; plotFaces(G, any([f1, f2], 2), 'edgecolor', 'r')
        
        % Schedule
        tsim = 1000*year;
        trep = [[1 5 10 30 60 180 365]*day [2 5 10:10:100 125:25:200 ...
                250:50:500 600:100:1000]*year];
        timesteps = [trep(1) diff(trep)];
        assert(sum(timesteps)== tsim, 'sum of timesteps must equal simTime')
        schedule = simpleSchedule(timesteps, 'W', [], 'bc', bc); 
        
        % Simulation
        N = 2;
        maxNumCompThreads(N);
        nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
        [~, states] = simulateScheduleAD(state0, model, schedule, ...
                                         'NonLinearSolver', nls, ...
                                          'Verbose', true);
        ds = norm(states{end-1}.s(:,1) - states{end}.s(:,1), inf);
        assert(ds < 1e-3)
        
        % Plot
        model = model.validateModel();
        for n=1:numel(states)
            states{n}.dp = states{n}.pressure - state0.pressure;
            states{n}.FlowProps.ComponentPhaseDensity = model.getProp(states{n}, 'CapillaryPressure'); % 2 components in brine phase.
            states{n}.FlowProps.CapillaryPressure = model.getProp(states{n}, 'CapillaryPressure');
            states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
        end
        states{1}.reg = model.rock.regions.saturation;
        
        figure(100)
        plotToolbar(G, states, 'edgealpha', 0.2)
        colormap(hot), colorbar;
    end
    
    % Create fcn
    sat_min=s_val(1);
    sat_max=s_val(end);
    pc_fcn=@(s) pcUpscaled(s, s_val, pc_val);
    
    % Output for table
    if sat_min == 0
        %sg = [sat_min 1e-3 0.01:0.01:sat_max];
        s_perc = s_val(find(s_val>s_pce, 1));
        sg = [sat_min 1e-5 linspace(s_perc, sat_max, opt.nval-2)];
    else
        error('check sat_min and adjust sg')
    end
    pc = pc_fcn(sg);
    
    
    % Plots?
    if nargin == 5 && makeplots == 1
        sg2 = [0 1e-3 1e-2:1e-2:sat_max];
        latx = {'Interpreter', 'latex'};
        figure(5)
        semilogy(sg, pc_fcn(sg)/1e5, '.-k')
        hold on
        if strcmp(opt.fault, 'test')
            colr = turbo(2);
            ylim([10^-3 10^3])
        else
            colr = turbo(numel(fluid.isclay));
            ylim([10^-1 10^3])
        end
        id_reg = unique(rock.regions.saturation);
        for n=1:size(colr, 1)
            semilogy(sg2, fluid.pcOG{id_reg(n)}(sg2)/1e5, '-', 'color', colr(n, :))
        end
        hold off
        grid on
        ylabel('$P_\mathrm{c}$ [bar]', 'fontsize', 14, latx{:})
        xlabel('$S_\mathrm{g}$ [-]', 'fontsize', 14, latx{:})
        xlim([0 1]); xticks([0:0.1:1])
    end

end

function varargout = pcUpscaled(s, s_val, pc_val)
  varargout{1}    = interpTable(s_val, pc_val, s);

  if nargout > 1
     varargout{2} = dinterpTable(s_val, pc_val, sg);
  end
end

function [model, states] = runCappaRutqvist2D(varargin)
    kf = 1e-16;
    mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

    % setup default option values
    opt = struct('method'             , 'fully coupled' , ...
                 'fluid_model'        , 'water', ...
                 'inj_rate'           , 0.02, ...
                 'mech_BC'            , 'left bottom roller', ...  
                 'nonlinearTolerance' , 1e-6            , ...
                 'splittingTolerance' , 1e-3            , ...
                 'verbose'            , false           , ...
                 'splittingVerbose'   , false           , ...
                 'properties'         , 'layers'        , ...
                 'grid'               , 'CR');
    opt = merge_options(opt, varargin{:});

    gravity reset y on;
    g = norm(gravity);

    %% 1. Mesh
    depth = 500; % depth of aquifer top surface
    Lx = 2000;Ly = 2000;
    if strcmp(opt.grid, 'uniform')
        nx = 50;
        ny = nx;
        G = cartGrid([nx, ny], [Lx Ly]); 
        G.nodes.coords(:,2) = G.nodes.coords(:,2) + depth;
        G = computeGeometry(G);
    elseif strcmp(opt.grid, 'CR')
        location_f = [500,1500];% [fault center x, fault center z]
        dip = 80;% fault dip (m)
        fW = 2.5;% fault width (m)
        f_D = 125; % fault displacement/throw window length (m)
        h_reservoir = 100;% reservoir thickness (m)
        h_caprock = 150;% caprock thickness (m)
        x_bottom = 500-cotd(dip)*Ly/2-fW/2;
        x_top = 500+cotd(dip)*Ly/2-fW/2;
        left_StretchFact_bottom = x_bottom/(location_f(1)-fW/2);
        left_StretchFact_top = x_top/(location_f(1)-fW/2);
        right_StretchFact_bottom = (Lx-x_bottom-fW)/(Lx-location_f(1)-fW/2);
        right_StretchFact_top = (Lx-x_top-fW)/(Lx-location_f(1)-fW/2);
            
        nxL = 20;nxR = 60;
        x = [linspace(0,location_f(1)-fW/2,nxL),linspace(location_f(1)+fW/2,Lx,nxR)];
        
        dy = fW*10; y = 500:dy:500+Ly;ny = length(y);
        G = tensorGrid(x,y);
        G = computeGeometry(G);
        
        fault_indx = find(G.cells.centroids(:,1)<location_f(1)+fW/2 & G.cells.centroids(:,1)>location_f(1)-fW/2);
        reservoir_indx = [find(G.cells.centroids(:,1)<location_f(1)-fW/2 &G.cells.centroids(:,2)>location_f(2)-h_reservoir/2 &G.cells.centroids(:,2)<location_f(2)+h_reservoir/2);...
                          find(G.cells.centroids(:,1)>location_f(1)+fW/2 &G.cells.centroids(:,2)>location_f(2)-h_reservoir/2-f_D &G.cells.centroids(:,2)<location_f(2)+h_reservoir/2-f_D)];
        caprock_indx = [find(G.cells.centroids(:,1)<location_f(1)-fW/2 &G.cells.centroids(:,2)>location_f(2)-h_reservoir/2-h_caprock &G.cells.centroids(:,2)<location_f(2)-h_reservoir/2);...
                        find(G.cells.centroids(:,1)<location_f(1)-fW/2 &G.cells.centroids(:,2)>location_f(2)+h_reservoir/2 &G.cells.centroids(:,2)<location_f(2)+h_reservoir/2+h_caprock);...
                        find(G.cells.centroids(:,1)>location_f(1)+fW/2 &G.cells.centroids(:,2)>location_f(2)-h_reservoir/2-h_caprock-f_D &G.cells.centroids(:,2)<location_f(2)-h_reservoir/2-f_D);...
                        find(G.cells.centroids(:,1)>location_f(1)+fW/2 &G.cells.centroids(:,2)>location_f(2)+h_reservoir/2-f_D &G.cells.centroids(:,2)<location_f(2)+h_reservoir/2+h_caprock-f_D)];
        basal_aquifer_indx = [find(G.cells.centroids(:,1)<location_f(1)-fW/2 &G.cells.centroids(:,2)>location_f(2)+h_reservoir/2+h_caprock);...
                          find(G.cells.centroids(:,1)>location_f(1)+fW/2 &G.cells.centroids(:,2)>location_f(2)+h_reservoir/2+h_caprock-f_D)];
                
        makeSkewL = @(c) ((left_StretchFact_top-left_StretchFact_bottom)/Ly*(Ly+depth-c(:,2))+left_StretchFact_bottom).*c(:,1);
        left_indx = find(G.nodes.coords(:,1)<= location_f(1)-fW/2);
        
        makeSkewR = @(c) Lx-((right_StretchFact_top-right_StretchFact_bottom)/Ly*(Ly+depth-c(:,2))+right_StretchFact_bottom).*(Lx-c(:,1));
        right_indx = find(G.nodes.coords(:,1)>= location_f(1)+fW/2);
        
            
        G.nodes.coords(left_indx,1) = makeSkewL(G.nodes.coords(left_indx,:));
        G.nodes.coords(right_indx,1) = makeSkewR(G.nodes.coords(right_indx,:));
        
    
    %     figure
    %     plotGrid(G,'FaceColor','none');
    %     plotGrid(G,fault_indx,'FaceColor','g');
    %     plotGrid(G,reservoir_indx,'FaceColor','b');
    %     plotGrid(G,caprock_indx,'FaceColor','y');
    %     set(gca, 'YDir','reverse');
    end

    %% 2. Rock Property
    % Here we define the porosity and permeability [L^2] of each grid cell. 
    poro = 0.1;
    rock.poro = repmat(poro, G.cells.num, 1);
    perm = 1e-14;
    rock.perm = repmat(perm, G.cells.num, 2);
    
    if strcmp(opt.properties, 'layers')
        assert(strcmp(opt.grid, 'CR'))
        rock.poro(caprock_indx) = 0.01;
    
        rock.perm(reservoir_indx,:) = repmat(1e-13*ones(1,2),size(reservoir_indx));
        rock.perm(caprock_indx,:) = repmat(1e-20*ones(1,2),size(caprock_indx));
        rock.perm(basal_aquifer_indx,:) = repmat(1e-15*ones(1,2),size(basal_aquifer_indx));
    
        if max(size(kf)) == 1
            rock.perm(fault_indx,:) = repmat([kf kf],size(fault_indx)); 
        else
            kf = reshape(kf,[2,5]);kf = kf';
            rock.perm(fault_top_indx,:) = repmat(1e-14*ones(1,2),size(fault_top_indx));
            rock.perm(fault_bot_indx,:) = repmat(1e-16*ones(1,2),size(fault_bot_indx));
            rock.perm(fault_throw1_indx,:) = repmat(kf(1,:),size(fault_throw1_indx));
            rock.perm(fault_throw2_indx,:) = repmat(kf(2,:),size(fault_throw2_indx));
            rock.perm(fault_throw3_indx,:) = repmat(kf(3,:),size(fault_throw3_indx));
            rock.perm(fault_throw4_indx,:) = repmat(kf(4,:),size(fault_throw4_indx));
            rock.perm(fault_throw5_indx,:) = repmat(kf(5,:),size(fault_throw5_indx));
        end
    end

%     figure
%     plotCellData(G,rock.poro,'EdgeColor','k');
%     colorbar;
%     set(gca, 'YDir','reverse');
%     
%     
%     figure
%     plotCellData(G,log10(rock.perm(:,2)),'EdgeColor','k');
%     colorbar;
%     set(gca, 'YDir','reverse');
% 
% 
    %% 3. Setup fluid parameters
    switch opt.fluid_model
      case 'blackoil-ccs'
        assert(strcmp(opt.grid, 'CR'))
        % Define 3 Saturation regions
        rock.regions.saturation = ones(G.cells.num, 1);        % upper aquifer
        rock.regions.saturation(caprock_indx) = 2;              % caprock
        rock.regions.saturation(fault_indx) = 3;                % fault
        rock.regions.rocknum = ones(G.cells.num,1); 
        fluid_path = fullfile(mrstPath('ls-proj'), '/uq-hannah/CappaRutqvist2D');
        fn  = fullfile(fluid_path, 'fPcFault_co2brine_modif_3ph_test.DATA');
        deck = convertDeckUnits(readEclipseDeck(fn));
        deck.REGIONS.ROCKNUM = rock.regions.rocknum;
        fluid = initDeckADIFluid(deck);
%         %% to be removed 
%         fluid.krW = fluid.krW{3};
%         fluid.krOW = fluid.krOW{3};
%         fluid.krG = fluid.krG{3};
%         fluid.krOG = fluid.krOG{3};
%         fluid.pcOG = fluid.pcOG{3};
      case 'blackoil-wg'
        pth = getDatasetPath('spe1');
        fn  = fullfile(pth, 'BENCH_SPE1.DATA');
        deck = readEclipseDeck(fn);
        deck = convertDeckUnits(deck);
        fluid = initDeckADIFluid(deck);
        if isfield(fluid, 'pcOW')
            fluid = rmfield(fluid, 'pcOW');
        end
        if isfield(fluid, 'pcOG')
            fluid = rmfield(fluid, 'pcOG');
        end
        % Setup quadratic relative permeabilities, since SPE1 relperm are a bit rough.
        fluid.krW = @(s) s.^2;
        fluid.krG = @(s) s.^2;
        fluid.krOW = @(s) s.^2;
        fluid.krOG = @(s) s.^2;
      case 'oil water'
        fluid = initSimpleADIFluid('phases', 'WO', 'mu', [1, 10]*centi*poise, ...
                                   'n',  [1, 1], 'rho', [1000, 700]*kilogram/ ...
                                   meter^3, 'c', 1e-10*[1, 1], 'cR', 4e-10);

      case 'water'
        fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
                                   1000, 'c', 1*1e-10, 'cR',4e-10);
      otherwise
        error('fluid_model  not recognized.');
    end

    %% 4. Setup material parameters for Biot and mechanics
    E          = 10 * giga * Pascal; % Young's module
    nu         = 0.25;               % Poisson's ratio
    alpha      = 1;                 % Biot's coefficient
    rock.alpha = alpha * ones(G.cells.num, 1);
    bulk_density = 2260*kilogram/meter^3;

    %% 5.1. ICs and BCs for flow
    f_top = find(G.faces.centroids(:,2) == min(G.faces.centroids(:,2)));
    f_bot = find(G.faces.centroids(:,2) == max(G.faces.centroids(:,2)));
    f_right = find(G.faces.centroids(:,1) == max(G.faces.centroids(:,1)));
    %f_left = find(G.faces.centroids(:,1) == min(G.faces.centroids(:,1)));
    
    f = [f_top;f_bot;f_right];
    [y_top, y_bot] = deal(min(G.faces.centroids(f,2)), max(G.faces.centroids(f,2)));

    clear initState;
    switch opt.fluid_model
      case 'blackoil-ccs'
        rho_or = fluid.rhoOS*kilogram/meter^3;
        p_top = 1*barsa + g*rho_or*y_top; % p at y_top
        equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [y_top, y_bot], p_top);
        fp_val = reshape(deval(equil, G.faces.centroids(f,2)), [], 1);
        bc = addBC([], f, 'pressure', fp_val, 'sat', [0, 1, 0]);
        initState.pressure = reshape(deval(equil, G.cells.centroids(:,2)), [], 1);
        initState.s  = repmat([0, 1, 0], [G.cells.num, 1]);
        initState.rs  = 0*fluid.rsSat(initState.pressure);
        clear equil;
      case 'blackoil-wg'
        rho_wr = fluid.rhoWS*kilogram/meter^3;
        p_top = 1*barsa + g*rho_wr*y_top; % p at y_top
        equil  = ode23(@(z,p) g *fluid.rhoWS(1), [y_top, y_bot], p_top);
        fp_val = reshape(deval(equil, G.faces.centroids(f,2)), [], 1);
        
        bc = addBC([], f, 'pressure', fp_val, 'sat', [1, 0, 0]);
        initState.pressure = reshape(deval(equil, G.cells.centroids(:,2)), [], 1);
        initState.s  = repmat([1, 0, 0], [G.cells.num, 1]);
        initState.rs  = 0*fluid.rsSat(initState.pressure);
        clear equil;
      case 'oil water'
        rho_or = fluid.rhoOS*kilogram/meter^3;
        p_top = 1*barsa + g*rho_or*y_top; % p at y_top
        equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [y_top, y_bot], p_top);
        fp_val = reshape(deval(equil, G.faces.centroids(f,2)), [], 1);
        bc = addBC([], f, 'pressure', fp_val, 'sat', [0, 1]);
        initState.pressure = reshape(deval(equil, G.cells.centroids(:,2)), [], 1);
        initState.s  = repmat([0, 1], [G.cells.num, 1]); 
        clear equil;
      case 'water'
        rho_wr = fluid.rhoWS*kilogram/meter^3;
        p_top = 1*barsa + g*rho_wr*y_top; % p at y_top
        equil  = ode23(@(z,p) g .* fluid.bW(p,0,false)*fluid.rhoWS(1), [y_top, y_bot], p_top);
        fp_val = reshape(deval(equil, G.faces.centroids(f,2)), [], 1);
        bc = addBC([], f, 'pressure', fp_val, 'sat', 1);
        initState.pressure = reshape(deval(equil, G.cells.centroids(:,2)), [], 1);
        initState.s  = ones(G.cells.num, 1); 
        clear equil;
      otherwise
        error('fluid_model not recognized.')
    end
  
%     figure 
%     plotCellData(G, initState.pressure/barsa, 'edgealpha', 0.2);
%     colormap(jet), c = colorbar;
%     c.Label.Interpreter = 'latex'; c.Label.FontSize = 11;
%     c.Label.String = '$p_0 $ [bar]';
%     set(gca, 'YDir','reverse');


    %% 5.2 BCs for mechanics
    switch opt.mech_BC
      case 'left bottom roller'
        bottom_nodes = find(G.nodes.coords(:,2) ==max(G.nodes.coords(:,2)));
        left_nodes = find(G.nodes.coords(:,1) ==min(G.nodes.coords(:,1)));
        
        left_bottom_nodes = intersect(left_nodes,bottom_nodes);
        left_nodes = setdiff(left_nodes,left_bottom_nodes);
        bottom_nodes = setdiff(bottom_nodes,left_bottom_nodes);
        
        disp_bc.nodes = [bottom_nodes;left_nodes;left_bottom_nodes];
        disp_bc.uu = repmat([0,0],G.nodes.num,1);
        disp_bc.mask = [repmat([false,true],numel(bottom_nodes),1);...% roller btm
                        repmat([true,false],numel(left_nodes),1);... % roller left
                        repmat([true,true],numel(left_bottom_nodes),1)];% left bottom corner fixed
        
        %% impose a given pressure at top and right, right sigma_h = 0.7*sigma_v
        pb_top = 1*barsa + g*bulk_density*y_top;
        equil  = ode23(@(z,p) g*bulk_density, [y_top, y_bot], pb_top); 
        pb_right_val = reshape(deval(equil, G.faces.centroids(f_right,2)), [], 1);
        clear equil;
        force = [repmat([0, pb_top],numel(f_top),1); 
                 [-0.7*pb_right_val,zeros(length(f_right),1)]];
        force_bc = struct('faces', [f_top;f_right], 'force', force);
    %% other cases not implemented yet    
    otherwise
        error('mech_BC not recognized.')
    end
    el_bc = struct('disp_bc' , disp_bc,'force_bc', force_bc);   

    %% Setup load for mechanics
    loadfun = @(x) repmat(bulk_density* [0 g], size(x, 1), 1);

    %% 6. Gather all the mechanical parameters in a struct
    % Transform these global properties (uniform) to cell values.
    E_matrix          = repmat(E, G.cells.num, 1);
    if strcmp(opt.grid, 'CR')
        E_matrix(fault_indx) = 5;
    end
    mech = struct('E', E_matrix, 'nu', repmat(nu, G.cells.num, 1), 'el_bc', el_bc, 'load', loadfun);

    %% 7. Setup model
    modeltype = [opt.method, ' and ', opt.fluid_model];
    fullycoupledOptions = {'verbose', opt.verbose};
    splittingOptions = {'splittingTolerance', opt.splittingTolerance, ...
                        'splittingVerbose', opt.splittingVerbose};
    switch modeltype
      case {'fully coupled and water'}
        model = MechWaterModel(G, rock, fluid, mech, fullycoupledOptions{:});
      case 'fully coupled and oil water'
        model = MechOilWaterModel(G, rock, fluid, mech, fullycoupledOptions{:});
      case {'fully coupled and blackoil-ccs','fully coupled and blackoil-wg'}
        model = MechBlackOilModel(G, rock, fluid, mech, fullycoupledOptions{:});
      case {'fixed stress splitting and blackoil-ccs','fixed stress splitting and blackoil-wg'}
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'blackoil', ...
                                               splittingOptions{:});
      case 'fixed stress splitting and oil water'
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'oil water', ...
                                               splittingOptions{:});
      case {'fixed stress splitting and water'}
        model = MechFluidFixedStressSplitModel(G, rock, fluid, mech, ...
                                               'fluidModelType', 'water', ...
                                               splittingOptions{:});
      otherwise
        error('modeltype not recognized.');
    end


    %% 8. Setup wells, well cell indices in 'global' grid
    wc_global = false(G.cartDims);
    center_iz = ceil(ny/2);
    wc_global(1, center_iz) = true; 
    wc = find(wc_global(G.cells.indexMap));

    switch opt.fluid_model
      case 'blackoil-wg'
            W = addWell([ ], G, rock, wc, 'Name', 'I1', 'Dir', 'z', ...
                        'Type', 'rate', 'Val', opt.inj_rate/fluid.rhoGS, 'compi', [0, 0, 1], ...    
                        'refDepth', G.cells.centroids(wc, G.griddim), ... 
                         'Radius', 0.2);
      case 'blackoil-ccs'
            W = addWell([ ], G, rock, wc, 'Name', 'I1', 'Dir', 'z', ...
                        'Type', 'rate', 'Val', opt.inj_rate/fluid.rhoGS, 'compi', [0, 0, 1], ...    
                        'refDepth', G.cells.centroids(wc, G.griddim), ... 
                         'Radius', 0.2);
      case 'oil water'
          W = addWell([ ], G, rock, wc, 'Name', 'I1', 'Dir', 'z', ...
                'Type', 'rate', 'Val', opt.inj_rate/fluid.rhoWS, 'compi', [1, 0], ...    
                'refDepth', G.cells.centroids(wc, G.griddim), ... 
                 'Radius', 0.2);
      case 'water'
            W = addWell([ ], G, rock, wc, 'Name', 'I1', 'Dir', 'z', ...
                'Type', 'rate', 'Val', opt.inj_rate/fluid.rhoWS, 'compi', 1, ...    
                'refDepth', G.cells.centroids(wc, G.griddim), ... 
                 'Radius', 0.2);
      otherwise
        error('fluid_model not recognized.')
    end

    facilityModel = FacilityModel(model.fluidModel);
    facilityModel = facilityModel.setupWells(W);
    model.FacilityModel = facilityModel;

    %% 9. Setup schedule
    reportTimes = [(12:12:60/minute)*minute, ...
           (2:1:24)*hour, ...
           (1440+5:5:1465)*minute, ... 
           ([26:2:48 54:6:96 108 120])*hour, ...
           [6:30 40:10:90]*day];
    timesteps = [reportTimes(1) diff(reportTimes)];
    assert(sum(timesteps)==reportTimes(end), 'sum of timesteps must equal simTime');
    schedule.control    = struct('W', W,'bc',bc);
    schedule.step.val = timesteps;
    schedule.step.control = [ones(numel(timesteps), 1)];
    
    %% compute mech IC
    % initial nodal displacements and pore pressure values are all zero
    % only need to determine the number of displacement values that must be
    % provided.
    initState.xd = zeros(nnz(~model.mechModel.operators.isdirdofs), 1);
    % add additional fields to initState, including u and uu (full lists of
    % displacement values, including the fixed displacement values associated
    % with imposed boundary conditions), as well as vdiv, stress and strain
    initState = addDerivedQuantities(model.mechModel, initState);
    

    %% 10. Simulate
    solver = NonLinearSolver('maxIterations', 100);
    [~, states, ~] = simulateScheduleAD(initState, model, schedule, 'nonlinearsolver', ...
                                        solver);
    model = model.validateModel();
    for n=1:numel(states)
        states{n}.dp = states{n}.pressure - initState.pressure;  
    %         pc = model.getProp(states{n}, 'CapillaryPressure');
    %         states{n}.FlowProps.CapillaryPressure = pc{1,2};
    %         states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, ...
    %                                                         'RelativePermeability');
    end
end
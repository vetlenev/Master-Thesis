function [s1D, vpar, G1D] = dynamicBLrun(G2, rock2, fluid, state0, ...
                            rate, ts, opt, states_plot)
        %
        %
        %
        solver_type = 'ad-bo';   % 'ad-bo' or 'incompr'
        
        % Grid and rock, same for both models 
        physDim = max(G2.nodes.coords) - min(G2.nodes.coords);
        G1D = cartGrid([G2.cartDims(end), 1], [physDim(end), physDim(1)]);
        %zmin = min(G2.nodes.coords(:,3));
        %G1D.nodes.coords(:,end) = G1D.nodes.coords(:,end) + zmin;
        G1D = computeGeometry(G1D);
        rock1D = makeRock(G1D, mean(rock2.perm(rock2.regions.saturation==1,:)), ...
                               mean(rock2.poro(rock2.regions.saturation==1)));
        gcNum = prod(G2.cartDims) - prod(G2.cartDims(1:2));
        rock1D.perm(1) = mean(rock2.perm(gcNum+1:end, end));
        rock1D.poro(1) = mean(rock2.poro(gcNum+1:end));
            
        if strcmp(solver_type, 'incompr')
            mrstModule add incomp
            p0m = mean(state0.pressure);
            mu = [fluid.muO(p0m) fluid.muG(p0m)];
            rho = [fluid.bO(p0m)*fluid.rhoOS fluid.bG(p0m)*fluid.rhoGS];
            fluid1D = initCoreyFluid('mu' , [mu(1),    mu(2)], ...
                'rho', [rho(1),  rho(2)], ...
                'n'  , [   5,      1.5], ...
                'sr' , [ 0.2,      0], ...
                'kwm', [   1,      1]);
            bc1D = fluxside([],  G1D, 'Left', rate/physDim(2), 'sat',  [0 1]);
            bc1D = pside(bc1D, G1D, 'Right', min(state0.pressure), 'sat', [1 0]);
            hT = computeTrans(G1D, rock1D);
            rSol = initState(G1D, [], unique(state0.pressure), [1 0]);
            rSol = incompTPFA(rSol, G1D, hT, fluid1D, 'bc', bc1D);

            % Tranport solver
            s1D = zeros(G1D.cells.num, numel(ts));
            rSole = cell(numel(ts), 1);
            for n=1:numel(ts)
                if n == 1
                    rSole{n} = implicitTransport(rSol, G1D, ts(n), rock1D, fluid1D, ...
                        'bc', bc1D, 'verbose', true);
                else
                    rSole{n} = implicitTransport(rSole{n-1}, G1D, ts(n), rock1D, ...
                        fluid1D, 'bc', bc1D, 'verbose', true);
                end
                s1D(:,n) = rSole{n}.s(:,2);
            end
            
        elseif strcmp(solver_type, 'ad-bo') 
            % Rock
            rock1D.regions.rocknum = ones(G1D.cells.num, 1);
            rock1D.regions.saturation = ones(G1D.cells.num, 1);
            
            % Fluid
            %       kw exp, kg exp,  Swc
            vpar = {3:.2:7, 1:.2:4, ...
                    fluid.krPts.og(1,2):0.02:fluid.krPts.og(2,2)};
            ntot = numel(vpar{1})*numel(vpar{2})*numel(vpar{3});
            kw_exp = repmat(vpar{1}', ntot/numel(vpar{1}), 1);
            kg_exp = repmat(repelem(vpar{2}, numel(vpar{1}))', ntot/(numel(vpar{1})*numel(vpar{2})), 1);
            swc = repelem(vpar{3}, numel(vpar{1})*numel(vpar{2}))';
            clear vpar;
            vpar = [kw_exp, kg_exp, swc];
            
            % Initial state
            uvp = unique(state0.pressure);
            uvp_d = diff(uvp)/mean(uvp);
            uvp(uvp_d < 1e-10) = [];
            state01d.pressure = uvp;
            state01d.s = repmat([1 0], G1D.cells.num, 1);
            state01d.rs = 0;
            state01d.rv = 0;
            
            s1D = cell(ntot, 1);
            %idp = 1:11:ntot;
            tic
            parfor n=1:ntot
                fluid1D = initSimpleADIFluid('phases', 'WG', 'n', [kw_exp(n) kg_exp(n)], ...
                                             'smin', [swc(n), 0], ...
                                             'mu',   [1 1], ...
                                             'rho',  [fluid.rhoOS fluid.rhoGS]);
                fluid1D.bG = fluid.bG;  fluid1D.bW = fluid.bO;
                fluid1D.muG = fluid.muG; fluid1D.muW = fluid.muO;
                fluid1D.pvMultR = fluid.pvMultR;

                % Model
                model = GenericBlackOilModel(G1D, rock1D, fluid1D, 'disgas', false, ...
                                             'vapoil', false, 'oil', false);

                % Acceleration and solver parameters (not needed)
                model = model.validateModel();

                % Model changes (after acceleration!)
                model.operators.p0 = state01d.pressure;                            % Needed for MyPvMult
                model.PVTPropertyFunctions = model.PVTPropertyFunctions.setStateFunction('PoreVolume', MyPvMult(model));
                model.minimumPressure = min(state01d.pressure);

                % BC - inject exact same equivalent volume as full 3D model
                bc1D = fluxside([],  G1D, 'Left', rate/physDim(2), 'sat',  [0 1]);
                bc1D = pside(bc1D, G1D, 'Right', min(state01d.pressure), 'sat', [1 0]);

                % Schedule
                tsim = opt.dyn_tsim_year*year;
                timesteps = [ts(1) diff(ts)];
                assert(sum(timesteps)== tsim, 'sum of timesteps must equal simTime')
                schedule = simpleSchedule(timesteps, 'W', [], 'bc', bc1D);

                % Simulate
                [~, states] = simulateScheduleAD(state01d, model, schedule, 'NonLinearSolver', ...
                                                 [], 'Verbose', false);                           
                s1D_it = cellfun(@(x) x.s(:,2), states, 'UniformOutput', false);
                s1D_it = reshape(cell2mat(s1D_it), G1D.cells.num, numel(ts));
                s1D_it = flipud(s1D_it);  % we want "top" cell first (zcor)
                s1D{n} = s1D_it;
                
                if mod(n, 100) == 0
                    disp('***********************************************')
                    disp([' Nsim = ' num2str(n) '/' num2str(ntot)]);
                    disp('***********************************************')
                end
            end
            toc
        end
        
        % Plot profiles
%         if states_plot
%             if strcmp(solver_type, 'incompr')
%                 states = rSole;
%             end
%             figure(randi(1000,1))
%             plotToolbar(G1D, states, 'edgealpha', 0.2); hold on
%             cmap = turbo;
%             colormap(cmap), c = colorbar; %clim([0 40000])
%             hold off
%             set(gca, 'ColorScale', 'linear')
%             caxis([0 1])
%             ax = gca;
%             ax.DataAspectRatio = [1 0.1 1];
%             %view([30 20])
%         end
        
    end
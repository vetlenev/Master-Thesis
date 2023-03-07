%% Perform CO2 injection (blackoil module)
%--------------------------------------------------------------------------
% medium FluidFlower analysis
%--------------------------------------------------------------------------

clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils 
mrstVerbose on

%% Options
temp_C     = 25; % TEMPERATURE in C (2525 for 25C_25std)
modelParams = 3; % 1 = initial, 2 = after area match, 3 = after full history match
inj_case   = 3;  % 1 = well test 1, 2 = well test 2, 3 = CO2inj
mesh_size  = 5;  % (fine) 1, 2, 4, 8 (coarse); refers to h (mm)
model_case = 3;  % 1 = medium model 1 props, 2 = medium model 2, 3 = medium model 3
inj_mult   = 1;  % injection rate multiplier
m = [];          % permeability multiplier (only for model_case 2)
removeS = true;  % remove cells corresponding to sealed part of silicone fault
dF = 3.0;        % diameter of F sand for perm computation (mm)
DiffCoeff = 1e-9;% diffusion coefficient [m2/s]
sgmin = 0.02;    % minimum saturation to activate relperm hysteresis
pcD = 0.3;       % Pc multiplier in unit D
revisedPc = 0;   % Pc revised after 1st run with experimental result at t=5h
flatMesh = false; % flatMesh (initial) or thickVar (depth-map based)
rampdown = 'fast';  % injection rampdown to 0. 'fast' (14s) or 'slow' (12min)
if model_case == 1
    folderName = ['T' num2str(temp_C) 'C_mParams' num2str(modelParams) '_bmesh' num2str(mesh_size) ...
                  '_thickVar' '_removeS' num2str(removeS) '_modelcase' num2str(model_case) ...
                 '_D' num2str(DiffCoeff) '_Im' num2str(inj_mult(1)) '_SgMinHyst' ...
                 num2str(sgmin) '_kdF' num2str(dF) 'mm' '_pcD' num2str(pcD)];
    if flatMesh
        folderName = ['exp' num2str(inj_case) '_bmesh' num2str(mesh_size) ...
                  '_flat' '_removeS' num2str(removeS) '_modelcase' num2str(model_case) ...
                  '_D' num2str(DiffCoeff) '_Im' num2str(inj_mult(1)) '_SgMinHyst' ...
                  num2str(sgmin) '_kdF' num2str(dF) 'mm' '_pcD' num2str(pcD) '_nohyst'];
    end
elseif model_case == 2
    folderName = ['T' num2str(temp_C) 'C_mParams' num2str(modelParams) '_bmesh' num2str(mesh_size) ...
                  '_thickVar' '_removeS' num2str(removeS) '_modelcase' num2str(model_case) ...
                  '_D' num2str(DiffCoeff) '_Im' num2str(inj_mult(1)) '_SgMinHyst' ...
                  num2str(sgmin)];
    if flatMesh
        folderName = ['exp' num2str(inj_case) '_bmesh' num2str(mesh_size) ...
                  '_flat' '_removeS' num2str(removeS) '_modelcase' num2str(model_case) ...
                  '_D' num2str(DiffCoeff) '_Im' num2str(inj_mult(1)) '_SgMinHyst' ...
                  num2str(sgmin)];
    end
elseif model_case == 3
    if revisedPc == 0
        folderName = ['T' num2str(temp_C) 'C_mParams' num2str(modelParams) '_bmesh' num2str(mesh_size) ...
                    '_thickVar' '_removeS' num2str(removeS) '_modelcase' num2str(model_case) ...
                    '_D' num2str(DiffCoeff) '_Im' num2str(inj_mult(1)) '_SgMinHyst' ...
                    num2str(sgmin)];
    elseif revisedPc == 1
        folderName = ['exp' num2str(inj_case) '_bmesh' num2str(mesh_size) ...
                      '_thickVar' '_removeS' num2str(removeS) '_modelcase' num2str(model_case) ...
                      '_D' num2str(DiffCoeff) '_Im' num2str(inj_mult(1)) '_SgMinHyst' ...
                      num2str(sgmin) 'revisedPcCDE'];
    elseif revisedPc == 2
        folderName = ['exp' num2str(inj_case) '_bmesh' num2str(mesh_size) ...
                      '_thickVar' '_removeS' num2str(removeS) '_modelcase' num2str(model_case) ...
                      '_D' num2str(DiffCoeff) '_Im' num2str(inj_mult(1)) '_SgMinHyst' ...
                      num2str(sgmin) 'revisedPcCE_nohyst'];
    end
    if flatMesh
        assert(revisedPc == 0)
        folderName = ['exp' num2str(inj_case) '_bmesh' num2str(mesh_size) ...
                     '_flat' '_removeS' num2str(removeS) '_modelcase' num2str(model_case) ...
                     '_D' num2str(DiffCoeff) '_Im' num2str(inj_mult(1)) '_SgMinHyst' ...
                     num2str(sgmin) '_nohyst'];
    end
end
%topDir = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\benchmark\';
if temp_C == 20
    topDir = '/home/lsalo/matlab/sim_data/mrst/fluidflower/benchmark/20C/';
elseif temp_C == 25
    topDir = '/home/lsalo/matlab/sim_data/mrst/fluidflower/benchmark/25C/';
elseif temp_C == 2525
    topDir = '/home/lsalo/matlab/sim_data/mrst/fluidflower/benchmark/25C_25std/';
end
%topDir = 'C:\Users\Lluis\matlab\sim_data\mrst\fluidflower\benchmark\';

[plts, opt] = bSimOpts(inj_case, mesh_size, inj_mult, removeS, ...
                       DiffCoeff, sgmin, dF, pcD, rampdown, flatMesh, ...
                       temp_C, modelParams);
if model_case == 2
    if opt.modelParams == 1
        m = [1 1 1 1 1 1 1];                % as provided by ex-situ measurements
    elseif opt.modelParams == 2
        m = [1 1 1 1 1 1 1];
    elseif opt.modelParams == 3             % matched model 2 perm multipliers
        m = [1 1 1 1 1.5 1.6 1];            % [ESF, ESFsup, C, D, E, F, G]          
    end
elseif model_case == 3
    if opt.modelParams == 1
        m = [1 1 1 1 1 1 1];
    elseif opt.modelParams == 2
        m = [1 1 1 1 1 1 1];
    elseif opt.modelParams == 3             % matched model 3
        if DiffCoeff == 1e-9
            m = [1/3 1/3 1 1 1.2 1.7 1.1]; 
        elseif DiffCoeff == 3e-9
            m = [1/3 1/3 1 1 1.6 2.3 1.3]; 
        end
    end
end


%% Mesh
[G, G_dat, unit, wellId, meshpt] = bMesh(opt, plts, removeS, 0);


%% Rock
rock = bRock(model_case, G, G_dat, unit, m, removeS, opt, plts, 0, 0);


%% Fluids
[fluid, rock] = bFluid(inj_case, model_case, G, G_dat, rock, unit, ...
                       removeS, revisedPc, opt, plts, 0);


%% Initialize
gravity reset on
[state0, p_r] = bInitialize(G, fluid, plts, 0);


%% Wells
[W, timesteps, wellInx] = bWells(inj_case, G, G_dat, rock, opt);


%% Model and acceleration
[model, nls] = bModel(G, rock, fluid, state0, inj_case, opt); 


%% BCs
[bc, model] = bBC(G, fluid, model, p_r, opt);


%% Schedule
schedule = bSchedule(timesteps, W, bc, inj_case, opt);


%% Simulation
N = 2;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                           'NonLinearSolver', nls, ...
%                                           'Verbose', true);
 
problem = packSimulationProblem(state0, model, schedule, folderName, 'Name', folderName, ...
                                'Directory', topDir, 'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem);  
% ^ Takes about 3 minutes to read data for benchmark 


%% Results
% Dense data
% Takes 10 minutes on GRS3
%ddata = getDenseData(timesteps, wellId, problem, unit, ...
%                     G_dat, states, wellSols, true);

% Areas map
% Takes about 3 minutes on GRS3
%[areas_A, areas_B] = getAreas(G, model, timesteps, states, problem, 1);

% Total Mass in domain (mobile, immobile, dissolved, total)
% Takes about 4 min on GRS3
%[total_mass] = getMass(model, states, timesteps, problem.Name, 1);

% Plots of Sg and Conc at given times
%       5  15  24  48  72  120h
tid = [314 442 496 640 784 1072];
%tid = 314;
sg_bound = 1e-3;
c_bound = 0.15*1.4;
plotsSgC(G, rock, model, meshpt, plts, states, sg_bound, c_bound, tid, timesteps)

% video
% dir = 'C:\Users\lsalo\matlab\';
% type = 'cgw';
% %idt = 1:314;
% %idt = [5:5:310 314 324:5:349 353:2:423 424:numel(timesteps)];
% idt = 1:numel(timesteps);
% %idt = [5:5:310 314 324:5:349 353:4:423 424:3:numel(timesteps)];
% bStatesVideo(dir, type, folderName, idt, plts, ...
%                G, rock, states, model, timesteps, meshpt);

% compute quantitites
model = model.validateModel();
bhp = zeros(numel(states), numel(W));
pw = zeros(numel(states), numel(wellId));
for n=1:numel(states)
    pw(n,:) = states{n}.pressure(wellId);
    states{n}.dp = states{n}.pressure - model.operators.p0;
    rho = model.getProp(states{n}, 'ComponentPhaseDensity');
    states{n}.FlowProps.ComponentPhaseDensity = rho{:,1}; % 2 components in brine phase.
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.FlowProps.CapillaryPressure = pc{1,2};
    states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
    bhp(n, :) = [wellSols{n}.bhp];
end
states{1}.reg = model.rock.regions.saturation;

% basic overview
%rescells = find(ismember(G_dat.p, unit.F));
%ESFcells = find(ismember(G_dat.p, unit.ESF));
plts.fig3D(); plotToolbar(G, states, 'edgealpha', 0.2); hold on
%plotGrid(G, rescells, 'edgecolor', [0.7 0.7 0.7]);
%plotGrid(G, ESFcells, 'edgecolor', [0.4 0.4 0.4]);
%plot(meshpt.boundary(:,1), meshpt.boundary(:,2), 'k')
xmx = max(G.faces.centroids(:,1));
for n=1:numel(meshpt.lines)
    if n < 5
        clr = [0.3 0.3 0.3];
    else
        clr = [0.7 0.7 0.7];
    end
    xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
    plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
end
%plotLinePath(meshpt.lines(5:end), 'b');
%plotLinePath(meshpt.lines(1:4), 'r');
cmap = flipud(cmocean('deep'));
plts.setAxProps(gca), colormap(cmap), c = colorbar; %clim([0 40000])
axis equal off
view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
set(gca, 'ColorScale', 'log')
caxis([1e-6 1])
ax = gca;
ax.DataAspectRatio = [0.1 1 1];

% Smooth pressure
tstep_injstop = find(cumsum(timesteps)/minute == 300) - 1;
pw2 = zeros(numel(timesteps), numel(wellId));
wl = 51;
for n=1:numel(timesteps)
    if n <= tstep_injstop
        id1 = max([1 n-(wl-1)/2]);
        id2 = min([n+(n-1) n+(wl-1)/2 tstep_injstop]);
        if n==1
            pw2(n,:) = pw(1,:);
        else
            if (tstep_injstop-n) < id1
                id1 = n - (tstep_injstop-n);
            end
            if n < tstep_injstop
                pw2(n,:) = mean(pw(id1:id2, :));
            else
                pw2(n,:) = pw(n, :);
            end
        end
        pw2(n, 2) = min([109560 pw2(n, 2)]);
        pw2(n, 3) = min([103700 pw2(n, 3)]);
    elseif n>tstep_injstop
        pw2(n,:) = model.operators.p0(wellId);
    end    
end

% Manual interp
tdat = [0 14 56 99 123 138 164 196 228 273 300 301 7200];
tid = [0 find(ismember(cumsum(timesteps)/minute, tdat))];
pw1_dat = zeros(numel(tdat), 1);
pw2_dat = zeros(numel(tdat), 1);
pw1_dat(1) = model.operators.p0(wellId(2));
pw2_dat(1) = model.operators.p0(wellId(3));
pw1_dat(2:end) = [109511 109518 109524 109528 109528 109531 109555 ...
                  109560 109494 109520 pw1_dat(1) pw1_dat(1)];        % pick values at half peak, then interp1
pw2_dat(2:end) = [103651 103660 103665 103667 103667 103671 103693 ...
                  103688 103641 103670 pw2_dat(1) pw2_dat(1)];
tdat(end) = 300;
tout = [0 10:10:7200];
pw1_out = interp1(tdat, pw1_dat, tout, 'pchip');
pw2_out = interp1(tdat, pw2_dat, tout, 'linear');

% Plot
subplot(1,2,1)
hold on; 
plot(cumsum(timesteps(1:300))/minute, pw(1:300, 2), '.-')
plot(cumsum(timesteps(1:312))/minute, pw2(1:312,2), 'linewidth', 1); 
plot(tout, pw1_out, 'linewidth', 1);
plot([tsm(142) tsm(142)], [109480 109600], '-k')
legend('Datapoints', 'Smooth', 'Interp', 'I2 start')
hold off
title('P1')
grid on
subplot(1,2,2)
hold on; 
plot(cumsum(timesteps(1:300))/minute, pw(1:300, 3), '.-')
plot(cumsum(timesteps(1:312))/minute, pw2(1:312,3), 'linewidth', 1); 
plot(tout, pw2_out, 'linewidth', 1);
plot([tsm(142) tsm(142)], [103620 103720], '-k')
hold off
title('P2')
grid on


% well cell pressure
dpw = zeros(numel(states), numel(wellId));
for n=1:numel(states)
    dpw(n, :) = states{n}.pressure(wellId) - model.operators.p0(wellId);
end
figure(2)
tsm = cumsum(timesteps(1:numel(states)))'/minute;
subplot(2,1,1)
hold on
plot(tsm(1:314), bhp(1:314,1)-model.operators.p0(wellId(1)), '-', 'color', [0.8 0.8 0.8], ...
     'linewidth', 2.5)
plot(tsm(143:314), bhp(143:314,2)-model.operators.p0(wellId(4)), '-', ...
     'color', [0.8 0.8 0.8], 'linewidth', 2.5)
plot(tsm, dpw)
text(305, 50, 'Injection stop')
legend('bh1', 'bh2', 'I1 cell', 'P1 cell', 'P2 cell', 'I2 cell')
xlabel('t [min]')
ylabel('\Delta p [Pa]')
grid on
xlim([0 sum(timesteps(1:numel(states)))/minute])
subplot(2,1,2)
hold on
plot(tsm(1:314), bhp(1:314,1)-model.operators.p0(wellId(1)), '-', 'color', [0.8 0.8 0.8], ...
     'linewidth', 2.5)
plot(tsm(143:314), bhp(143:314,2)-model.operators.p0(wellId(4)), '-', ...
     'color', [0.8 0.8 0.8], 'linewidth', 2.5)
plot(tsm(1:314), dpw(1:314, :))
%legend('bh1', 'bh2', 'I1 cell', 'P1 cell', 'P2 cell', 'I2 cell')
xlabel('t [min]')
ylabel('\Delta p [Pa]')
grid on
xlim([0 sum(timesteps(1:314))/minute])



% CO2 mass and V end of I1 (not accounting for rock porevolume changes)
% tid = 54;
% sg_bound = 1e-2;
% V_with_co2g = sum(G.cells.volumes(states{tid}.s(:,2) >= sg_bound)); % Sum of cell volumes
% V_with_co2b = sum(G.cells.volumes(states{tid}.PVTProps.Density{1} > ...
%                           995.18));
% Vco2_g_cells = model.operators.pv .* states{tid}.s(:,2);
% Vco2_g = sum(Vco2_g_cells);                                     % Actual g volume
% Mco2_g = sum(Vco2_g_cells.*states{tid}.PVTProps.Density{2});
% Mco2_tot = sum(states{tid}.FlowProps.ComponentTotalMass{2});
% Mco2_b = Mco2_tot - Mco2_g;
% rsSat = fluid.rsSat(states{tid}.pressure);
% flag = states{tid}.rs >= rsSat;
% bO = fluid.bO(states{tid}.pressure, states{tid}.rs, flag);
% Vco2_b = sum(model.operators.pv .* states{tid}.s(:,1) .* states{tid}.rs .* bO);
% Vco2_tot = Vco2_b + Vco2_g;
% %Mco2_b_check = Vco2_b*fluid.rhoGS;
% disp(['Surface V of dissolved CO2 is ' num2str(Vco2_b*10^6) ' ml'])
% disp(['Mass of CO2 in the brine is ' num2str(Mco2_b*10^6) ' mg'])
% disp(['total V of CO2 in the domain is ' num2str(Vco2_tot*10^6) ' ml'])
% disp(['total mass of CO2 in the domain is ' num2str(Mco2_tot*10^6) ' mg'])
% 
% % End of I2
% tid = 122;
% sg_bound = 1e-2;
% V_with_co2g2 = sum(G.cells.volumes(states{tid}.s(G_dat.p==4,2) >= sg_bound)); % Sum of cell volumes
% V_with_co2b2 = sum(G.cells.volumes(all([states{tid}.PVTProps.Density{1} > ...
%                                        995.15, G_dat.p==4], 2)));
%                                    
% V_with_co2g_mid = sum(G.cells.volumes(states{tid}.s(G_dat.p==2,2) >= 1e-3)); % Sum of cell volumes
% V_with_co2b_mid = sum(G.cells.volumes(all([states{tid}.PVTProps.Density{1} > ...
%                                        995.15, G_dat.p==2], 2)));

% Plot cell data at a given time
% idt = 88;
% plts.fig3D();
% plotCellData(G, states{idt}.s(:,2), 'edgecolor', 'none')
% %plotCellData(G, states{idt}.PVTProps.Density{1}, 'edgecolor', 'none')
% hold on
% xmx = max(G.faces.centroids(:,1));
% for n=1:numel(meshpt.lines)
%     if n < 5
%         clr = [0.3 0.3 0.3];
%     else
%         clr = [0.7 0.7 0.7];
%     end
%     xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
%     plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
% end
% plts.setAxProps(gca)
% colormap(hot)
% %colormap(cmocean('balance'))
% c = colorbar; 
% ylabel(c, '$S_{g}$ [-]', 'fontsize', 14, 'interpreter', 'latex')
% %ylabel(c, '$\rho_\mathrm{aq}$ [kg/m$^3$]', 'fontSize', 14, 'interpreter', 'latex')
% %caxis([995.1 995.4])
% axis equal off
% view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
% set(gca, 'ColorScale', 'log')
% caxis([1e-4 1])

% Brine Density at a given time
% id = states{end}.rs > 0.001;
% plotGrid(G,'edgealpha', 0.2)
% plotCellData(G, states{end}.PVTProps.Density{1}(id), find(id), 'edgealpha', 0.2)
% setAxProps(gca), colormap(turbo), c = colorbar;
% axis equal off, view([90 0])
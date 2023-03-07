%% Perform CO2 injection (blackoil module)
%--------------------------------------------------------------------------
% medium FluidFlower analysis
%--------------------------------------------------------------------------

clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils 
mrstVerbose off


%% Options
inj_case = 2;       % 1 = well test 1, 2 = well test 2
mesh_size =  4;     % (fine) 1, 2, 4, 8 (coarse); refers to h (mm)
model_case = 2;     % 1 = medium model 1 props, 2 = benchmark reported props
batch = 'second';    
removeS = true;
%folderName = ['exp' num2str(inj_case) '_bmesh' num2str(mesh_size) ...
%              '_modelcase' num2str(model_case) '_D1e-9' ...
%              '_mediumProps_FPc0_ESF0.1mm_ESFPce1e-6'];
%topDir = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\benchmark\';
topDir = 'C:\Users\Lluis\matlab\sim_data\mrst\fluidflower\benchmark\welltest\wt1_2\';

[plts, opt] = bSimOpts(inj_case, mesh_size, [], removeS);
m = [0.5, 1.5, 1.4, 1.3, 0.7, 0.7, 1; ...      % min MAE
     0.834, 1.5, 0.6, 0.7, 0.7, 0.7, 0.7];     % min STD
nsim = size(m, 1);

%% Mesh
[G, G_dat, unit, wellId, meshpt] = bMesh(opt, plts, removeS, 0);

for k=1:nsim
    
    folderName = ['exp' num2str(inj_case) '_bmesh' num2str(mesh_size) ...
        '_modelcase' num2str(model_case) ...
        '_reportedProps_stepMaxRate_kmult_ESF', num2str(m(k,1)) ...
        '_ESFsup', num2str(m(k,2)), '_C', num2str(m(k,3)), '_D', num2str(m(k,4)) ...
        '_E', num2str(m(k,5)), '_F', num2str(m(k,6)), '_G', num2str(m(k,7))];
    

    %% Rock
    rock = bRock(model_case, G, G_dat, unit, m(k, :), removeS, plts, 0, 0);
    
    
    %% Fluids
    [fluid, rock] = bFluid(inj_case, model_case, G, G_dat, rock, ...
                           unit, removeS, opt, plts, 0);
    
    
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
    N = 4;
    maxNumCompThreads(N);
    nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
    %[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
    %                                           'NonLinearSolver', nls, ...
    %                                           'Verbose', true);
    
    problem = packSimulationProblem(state0, model, schedule, folderName, 'Name', folderName, ...
        'Directory', topDir, 'NonLinearSolver', nls);
    disp('***************************************************************')
    disp(['Running batch ' num2str(batch)])
    disp('***************************************************************')
    [ok, status] = simulatePackedProblem(problem);
    [wellSol, states, report] = getPackedSimulatorOutput(problem);

end

%% Results
% compute quantitites
%model = model.validateModel();
bhp = zeros(numel(states), numel(W));
fports = [0.0120 0.53 1.01; 0.0120 0.53 0.61];
dist = pdist2(G.cells.centroids, fports);
[~, fpInx] = min(dist,[],1);
wellId = [wellId; fpInx'];
wellId_dp = zeros(numel(states), numel(wellId));
wellId_dbhp = zeros(numel(states), numel(W));
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;
    rho = model.getProp(states{n}, 'ComponentPhaseDensity');
    %states{n}.FlowProps.ComponentPhaseDensity = rho{:,1}; % 2 components in brine phase.
    %pc = models{1}.model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.FlowProps.CapillaryPressure = pc{1,2};
    %states{n}.FlowProps.RelativePermeability = models{1}.model.getProp(states{n}, 'RelativePermeability');
    bhp(n, :) = [wellSol{n}.bhp];
    wellId_dp(n, :) = (states{n}.pressure(wellId) - state0.pressure(wellId))';
    wellId_dbhp(n, :) = [wellSol{n}.bhp]' - state0.pressure([W.cells]);
end
states{1}.reg = model.rock.regions.saturation;


% basic overview
%rescells = find(ismember(G_dat.p, unit.F));
%ESFcells = find(ismember(G_dat.p, unit.ESF));
plts.fig3D(); plotToolbar(G, states, rock.regions.saturation ~= 7, 'edgealpha', 0); hold on
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
plts.setAxProps(gca), colormap(hot), c = colorbar; %clim([0 40000])
axis equal off
view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
% set(gca, 'ColorScale', 'log')
% caxis([1e-4 1])

% Load Well test data and compare with pressure at corresponding cell
pth = 'C:\Users\Lluis\matlab\mrst-dev\mrst-ls\ls-proj\fluidflower\benchmark\welltests_data\';
np = 9;
if inj_case == 1
    fn = 'pressure_measurements_welltest_port_17_7.xls';
    id_tcol = [7 7+(3:3:3*(np-1))];
    id_pcol = [8 8+(3:3:3*(np-1))];
    id_t0 = [repelem(1143, 6) 1142, 1143, 1142];
elseif inj_case == 2
    fn = 'pressure_measurements_welltest_port_9_3.xls';
    id_tcol = [4 4+(3:3:3*np)]; id_tcol(id_tcol==28) = [];
    id_pcol = [5 5+(3:3:3*np)]; id_pcol(id_pcol==29) = [];
    id_t0 = repelem(228, 9);
end
data = readtable([pth, fn], 'Sheet', 1);
p_data = cell(1, np);
t_data = cell(1, np);
for n=1:np
    tval = table2array(data(:, id_tcol(n)));
    t0 = tval(id_t0(n));
    tval = tval - t0;
    tfin = opt.schedule(end)*minute;
    [~, id_tfin] = min(abs(tfin - tval));
    t_data{n} = tval(id_t0(n):id_tfin);
    p_data{n} = table2array(data(id_t0(n):id_tfin, id_pcol(n)))*1e3;    % mbar
end
if inj_case == 1
    id_i1 = 4;
    id_i2 = 1;
    id_p1 = [5 7];
    id_p2 = [8 9];
    id_f1 = [2 6];
    id_f2 = 3;
elseif inj_case == 2
    id_i1 = 4;
    id_i2 = 8;
    id_p1 = [1 9];
    id_p2 = [5 7];
    id_f1 = [2 6];
    id_f2 = 3;
end
p_i2 = p_data{id_i2};         % 0-2.5 bar
p_i1 = p_data{id_i1};         % -1-2.5 bar
p_p1_1 = p_data{id_p1(1)};       % 0-6 bar
p_p1_2 = p_data{id_p1(2)};       % 0-4 bar
p_p2_1 = p_data{id_p2(1)};       % 0-100,
p_p2_2 = p_data{id_p2(2)};       % 0-4 bar
p_f1_1 = p_data{id_f1(1)};       % 0-2.5,
p_f1_2 = p_data{id_f1(2)};       % 0-100 bar
p_f2 = p_data{id_f2};         % 0-2.5 bar

% Plot
latx={'Interpreter','latex'};
prpl = [92, 16, 199]/255;
figure(37)
subplot(2,3,1)
    hold on
    p2 = plot(t_data{id_i1}/minute, p_i1, '--k', 'linewidth', 0.5, ...
              'DisplayName', '[-1--2.5] bar');
%     if inj_case == 1
        p1 = plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dp(:,1)/100, '-r', 'linewidth', 1, ...
                     'DisplayName', 'Simulation');
%     elseif inj_case == 2
%         p1 = plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dbhp/100, '-r', 'linewidth', 1, ...
%                      'DisplayName', 'Simulation');
%     end
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('I$_1$ $\vert$ $y, z = 0.93, 1.01$m $\vert$ [9, 3]', latx{:})
    grid on
    xlabel('$t$ [min]', latx{:})
    ylabel('$\Delta p$ [mbar]', latx{:})
    h=legend([p1, p2], 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,2)
    hold on
    p1 = plot(t_data{id_i2}/minute, p_i2, '--k', 'linewidth', 0.5, ...
         'DisplayName', '[0--2.5] bar');
%     if inj_case == 2
        plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dp(:,4)/100, '-r', 'linewidth', 1)
%     elseif inj_case == 1
%         plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dbhp/100, '-r', 'linewidth', 1)
%     end
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('I$_2$ $\vert$ $y, z = 1.73, 0.61$m $\vert$ [17, 7]', latx{:})
    grid on
    h=legend(p1, 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,3)
    hold on
    p1 = plot(t_data{id_p1(1)}/minute, p_p1_1, '--k', 'linewidth', 0.5, ...
         'DisplayName', '[0--6] bar');
    p2 = plot(t_data{id_p1(2)}/minute, p_p1_2, '--', 'color', [0.6 0.6 0.6], ...
         'linewidth', 0.5, 'DisplayName', '[0--4] bar');
    plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dp(:,2)/100, '-', 'color', prpl, ...
        'linewidth', 1)
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('P$_1$ $\vert$ $y, z = 1.53, 0.81$m $\vert$ [15, 5]', latx{:})
    grid on
    h=legend([p1, p2], 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,4)
    hold on
    p1 = plot(t_data{id_p2(1)}/minute, p_p2_1, '--k', 'linewidth', 0.5, ...
         'DisplayName', '[0--100] bar');
    p2 = plot(t_data{id_p2(2)}/minute, p_p2_2, '--', 'color', [0.6 0.6 0.6], ...
         'linewidth', 0.5, 'DisplayName', '[0--4] bar');
    plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dp(:,3)/100, '-', 'color', prpl, ...
         'linewidth', 1)
    grid on
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('P$_2$ $\vert$ $y, z = 1.73, 0.21$m $\vert$ [17, 11]', latx{:})
    hold off
    h=legend([p1, p2], 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,5)
    hold on
    p1 = plot(t_data{id_f1(1)}/minute, p_f1_1, '--k', 'linewidth', 0.5, ...
         'DisplayName', '[0--2.5] bar');
    p2 = plot(t_data{id_f1(2)}/minute, p_f1_2, '--', 'color', [0.6 0.6 0.6], ...
         'linewidth', 0.5, 'DisplayName', '[0--100] bar');
    plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dp(:,5)/100, '-c',  ...
         'linewidth', 1)
    grid on
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('F$_1$ $\vert$ $y, z = 0.53, 1.01$m $\vert$ [5, 3]', latx{:})
    hold off
    h=legend([p1, p2], 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,6)
    hold on
    p1 = plot(t_data{id_f2}/minute, p_f2, '--k', 'linewidth', 0.5, ...
         'DisplayName', '[0--2.5] bar');
    plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dp(:,6)/100, '-c', ...
         'linewidth', 1)
    grid on
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('F$_2$ $\vert$ $y, z = 0.53, 0.61$m $\vert$ [5, 7]', latx{:})
    hold off
    h=legend(p1, 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));

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
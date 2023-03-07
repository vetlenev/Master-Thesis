%% Convergence analysis
%
%
%
clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils 

% Load data (mesh, rock, fluid, states)
inj_type = 1;
model_case = 1;
mesh_sizes = [1 2 4 8];
ngrid = numel(mesh_sizes);
mesh_name = 'G_composite_cellSize';
problem_name = {'exp1_mesh1_model1_problem.mat', ...
                'exp1_mesh2_model1_problem.mat', ...
                'exp1_mesh4_model1_problem.mat', ...
                'exp1_mesh8_model1_problem.mat'};
sim_data_topfolder = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\medium\';
sim_data_folder = {'exp1_mesh1_modelcase1', ...
                   'exp1_mesh2_modelcase1', ...
                   'exp1_mesh4_modelcase1', ...
                   'exp1_mesh8_modelcase1'};
unit.ESF = [1 16 7];
unit.C = [3 5 11 13 15];
unit.E = [6 8 10 12];
unit.F = [2 4 9 14 17];
G_all = cell(ngrid,1);
states_all = cell(ngrid, 1);
rock_all = cell(ngrid, 1);
p_all = cell(ngrid, 1);
for n=1:4 %1:numel(mesh_size)
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/mesh');
    G_dat  = load(fullfile(pth, [mesh_name num2str(mesh_sizes(n)) '.mat']));  
    p_all{n} = G_dat.p;
    G_all{n} = G_dat.G;
    load([sim_data_topfolder problem_name{n}]);
    rock_all{n} = rock;
    [~, states_all{n}, ~] = getPackedSimulatorOutput(problem);  
    for k=1:numel(states_all{n})
        states_all{n}{k}.FlowProps.ComponentPhaseMass = ...
                    model.getProp(states_all{n}{k}, 'ComponentPhaseMass');
    end
end

% Clean extra tsteps from states in mesh 1

% Plotting
fig3D = @() figure('Position', [0, 0, 1300, 650]);
setAxProps = @(ax) set(ax, 'View'              , [65, 20]         , ...
                       'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                       'Projection'        , 'Perspective'     , ...
                       'Box'               , 'on');
latx={'interpreter', 'latex'};
                   

%% Free-phase CO2 mass vs time 
% (avg in inj reservoir 1 and 2)
% Avg saturation and plume length at top of each reservoir (y)
unit.inj1 = 9;
unit.inj2 = 4;
tsm = (cumsum(timesteps)/minute)';
tsh = (cumsum(timesteps)/hour)';
s_avg1 = zeros(numel(states_all{1}), ngrid);
s_avg2 = zeros(numel(states_all{1}), ngrid);
id_inj1{n} = cell(ngrid, 1);
id_inj2{n} = cell(ngrid, 1);
for n=1:4
    id_inj1{n} = ismember(p_all{n}, unit.inj1);
    id_inj2{n} = ismember(p_all{n}, unit.inj2);
    % fig3D(); plotGrid(G_all{3}, 'faceColor', 'none', 'edgealpha', 0.2);
    % plotGrid(G_all{3}, id_inj1, 'faceColor', 'none', 'edgecolor', 'r');
    % setAxProps(gca), axis equal; view([90 0])
  for k=1:numel(states_all{n})
      s_avg1(k, n) = mean(states_all{n}{k}.s(id_inj1{n}, 2));
      s_avg2(k, n) = mean(states_all{n}{k}.s(id_inj2{n}, 2));
  end
end

% Actual figure
fh = figure(1);
subplot(1,2,1)
hold on
plot(tsh(1:21), s_avg1(1:21,1), '-', 'color', [0.7 0.7 0.7], ...
     'linewidth', 3, 'DisplayName', '$h=1$ [mm]')
plot(tsh, s_avg1(:,2), '-k', 'linewidth', 1.5, 'DisplayName', '$h=2$')
plot(tsh, s_avg1(:,3), '-r', 'linewidth', 1, 'DisplayName', '$h=4$')
plot(tsh, s_avg1(:,4), '-b', 'linewidth', 0.5, 'DisplayName', '$h=8$')
xlabel('$t$ [h]', latx{:}, 'fontsize', 14)
ylabel('$\overline{S_g}$', latx{:}, 'fontSize', 14)
xlim([0 10]), xticks(0:10)
%set(gca,'xscale','log')
legend(latx{:}, 'fontsize', 11)
grid on
hold off

subplot(1,2,2)
hold on
plot(tsh(1:21), s_avg2(1:21,1), '-', 'color', [0.7 0.7 0.7], ...
     'linewidth', 3, 'DisplayName', '$h=1$ [mm]')
plot(tsh, s_avg2(:,2), '-k', 'linewidth', 1.5, 'DisplayName', '$h=2$')
plot(tsh, s_avg2(:,3), '-r', 'linewidth', 1, 'DisplayName', '$h=4$')
plot(tsh, s_avg2(:,4), '-b', 'linewidth', 0.5, 'DisplayName', '$h=8$')
xlabel('$t$ [h]', latx{:}, 'fontsize', 14)
ylabel('$\overline{S_g}$', latx{:}, 'fontSize', 14)
xlim([0 10]), xticks(0:20)
%set(gca,'xscale','log')
legend(latx{:}, 'fontsize', 11)
grid on
hold off
set(fh, 'position', [200, 200, 600, 300]);

% Fig 2 - plume length at top of reservoir (end of inj 1 and 2)
% Find missing cells (first run loop below to check if missing cells)
% plotGrid(G_all{n}); hold on; plotFaces(G_all{n}, f1, 'edgecolor', 'r')
% view([90 0])
% idp1 = G_all{n}.cells.centroids(:, 2) > 0.5;
% idp2 = G_all{n}.cells.centroids(:, 2) < 0.6;
% idp = all([idp1, idp2], 2);
% idp = any([all([id_e1 idp], 2) all([id_inj1{n} idp], 2)], 2);
% text(0.0105*ones(sum(idp),1), G_all{n}.cells.centroids(idp,2), ...
%      G_all{n}.cells.centroids(idp,3), num2str(find(idp)));
% ylim([0.5 0.6]), zlim([0.3 0.4])
% n=1, f1: 1149, 1165
% n=2, f1: 1040
% n=3, f1: 532 548, f2: []
% n=4,

% Find cells at top of each reservoir F
c_inj1_top=cell(ngrid, 1);
c_inj2_top=cell(ngrid, 1);
l = zeros(1, ngrid);
for n=1:4
    % top cells in inj1 F compartment 
    id_e1 = ismember(p_all{n}, unit.ESF(1));
    fe1 = any([ismember(G_all{n}.faces.neighbors(:,1), find(id_e1)) ...
               ismember(G_all{n}.faces.neighbors(:,2), find(id_e1))], 2);
    ff1 = any([ismember(G_all{n}.faces.neighbors(:,1), find(id_inj1{n})) ...
               ismember(G_all{n}.faces.neighbors(:,2), find(id_inj1{n}))], 2);
    f1 = all([fe1, ff1], 2); % dividing faces 
    c1 = unique(G_all{n}.faces.neighbors(f1,:));
    c_inj1 = find(id_inj1{n});
    id_top = ismember(c_inj1, c1);
    c_inj1_top{n} = c_inj1(id_top);
    if n == 1
        toadd = [1149 1165];
    elseif n == 2       % cells missing in faces.neighbors
        toadd = 1040;
    elseif n == 3
        toadd = [532 548];
    elseif n ==4
        toadd = [1149 1165];
    end
    c_inj1_top{n} = unique([toadd'; c_inj1(id_top)]);
    
    % top cells in inj2 F compartment 
    id_e2 = ismember(p_all{n}, unit.ESF(3));
    fe2 = any([ismember(G_all{n}.faces.neighbors(:,1), find(id_e2)) ...
               ismember(G_all{n}.faces.neighbors(:,2), find(id_e2))], 2);
    ff2 = any([ismember(G_all{n}.faces.neighbors(:,1), find(id_inj2{n})) ...
               ismember(G_all{n}.faces.neighbors(:,2), find(id_inj2{n}))], 2);
    f2 = all([fe2, ff2], 2); % dividing faces 
    c2 = unique(G_all{n}.faces.neighbors(f2,:));
    c_inj2 = find(id_inj2{n});
    id_top = ismember(c_inj2, c2);
    c_inj2_top{n} = c_inj2(id_top);       
    
    % plume length, tstep=18 (last at max rate before rampdown)
    id_g1 = states_all{n}{18}.s(c_inj1_top{n},2) > 1e-4;
    cc1 = G_all{n}.cells.centroids(c_inj1_top{n}, 2);
    l(n) = max(cc1(id_g1))-min(cc1(id_g1));
end

fh = figure(2);
hold on
plot(1./mesh_sizes, l*100, '-sk', 'markersize', 8, 'markerfacecolor', [0.5 0.5 0.5])
xlabel('$1/h$ [1/mm]', latx{:}, 'fontsize', 14)
ylabel('$L$ [mm]', latx{:}, 'fontSize', 14)
xlim([0 1.25]), xticks(0:0.25:1.25)
ylim([40 50])
%set(gca,'xscale','log')
grid on
hold off
set(fh, 'position', [200, 200, 350, 250]);


%% Residual (trapped) CO2 mass vs time
% compute sgt at each timestep, in cells of interest (inj 1 and inj 2)
% adjust sgt if dissolved
% mass = sgt*(poro*cell_vol)*co2dens [kg]
a = fluid.ehystr{4};        % from fluid.ehystr{4} in .DATA file
sgmx = fluid.krPts.g(5, 3); % max CO2 sat (F unit)
sgmn = 0;                   % no residual sat during primary drainage
sgtmx = fluid.sgtmax(1);
A = @(sgiv) 1 + a*(sgmx - sgiv);    % sgiv = sgmax in run
C = 1/(sgtmx - sgmn) - 1/(sgmx-sgmn);
sgt = @(sgiv) sgmn + (sgiv-sgmn)./(A(sgiv) + C*(sgiv-sgmn));
trapCO2m_i1 = zeros(numel(states_all{end}),ngrid);
trapCO2m_i2 = zeros(numel(states_all{end}),ngrid);
st_avg1 = zeros(numel(states_all{end}),ngrid);
st_avg2 = zeros(numel(states_all{end}),ngrid);
for n=2:4
   for k=1:numel(states_all{n})
       smax1 = states_all{n}{k}.sMax(:, 2);
       s1 = states_all{n}{k}.s(:, 2);
       id_sgmx = smax1 + 1e-3 >= sgmx;
       id_killo = ~id_sgmx;
       sgtv = zeros(numel(s1), 1);
       if any(id_sgmx)
           sgtv(id_sgmx) = fluid.sgtmax(1);
       end
       sgtv(id_killo) = sgt(smax1(id_killo));
       if any(sgtv > s1)
          sgtv(sgtv>s1) = s1(sgtv>s1); 
       end
       sgtv_i1 = sgtv(id_inj1{n});
       sgtv_i2 = sgtv(id_inj2{n});
       st_avg1(k,n) = mean(sgtv_i1);
       st_avg2(k,n) = mean(sgtv_i2);
       trapCO2m_i1(k,n) = sum(sgtv_i1.*rock_all{n}.poro(id_inj1{n}).*...
                              G_all{n}.cells.volumes(id_inj1{n}).*...
                              states_all{n}{k}.PVTProps.Density{2}(id_inj1{n}));
       trapCO2m_i2(k,n) = sum(sgtv_i2.*rock_all{n}.poro(id_inj2{n}).*...
                              G_all{n}.cells.volumes(id_inj2{n}).*...
                              states_all{n}{k}.PVTProps.Density{2}(id_inj2{n}));
   end
    
end

% Plot
fh = figure(3);
subplot(1,2,1)
hold on
%plot(tsh, st_avg1(:,1), '-', 'color', [0.7 0.7 0.7], ...
%     'linewidth', 2, 'DisplayName', '$h=1$ [mm]')
plot(tsh, st_avg1(:,2), '-k', 'linewidth', 1.5, 'DisplayName', '$h=2$')
plot(tsh, st_avg1(:,3), '-r', 'linewidth', 1, 'DisplayName', '$h=4$')
plot(tsh, st_avg1(:,4), '-b', 'linewidth', 0.5, 'DisplayName', '$h=8$')
xlabel('$t$ [h]', latx{:}, 'fontsize', 14)
ylabel('$\overline{S_{gt}}$', latx{:}, 'fontSize', 14)
xlim([0 10]), xticks(0:10)
%set(gca,'xscale','log')
legend(latx{:}, 'fontsize', 11)
grid on
hold off

subplot(1,2,2)
hold on
%plot(tsh, st_avg2(:,1), '-', 'color', [0.7 0.7 0.7], ...
%     'linewidth', 2, 'DisplayName', '$h=1$ [mm]')
plot(tsh, st_avg2(:,2), '-k', 'linewidth', 1.5, 'DisplayName', '$h=2$')
plot(tsh, st_avg2(:,3), '-r', 'linewidth', 1, 'DisplayName', '$h=4$')
plot(tsh, st_avg2(:,4), '-b', 'linewidth', 0.5, 'DisplayName', '$h=8$')
xlabel('$t$ [h]', latx{:}, 'fontsize', 14)
ylabel('$\overline{S_{gt}}$', latx{:}, 'fontSize', 14)
xlim([0 10]), xticks(0:20)
%set(gca,'xscale','log')
legend(latx{:}, 'fontsize', 11)
grid on
hold off
set(fh, 'position', [200, 200, 600, 300]);


%% Dissolved CO2 mass vs time
disCO2m_i1 = zeros(numel(states_all{end}), ngrid);
disCO2m_i2 = zeros(numel(states_all{end}), ngrid);
for n=2:4
    for k = 1:numel(states_all{n})
        disCO2m_i1(k, n) = sum(states_all{n}{k}.FlowProps.ComponentPhaseMass{2,1}(id_inj1{n}));
        disCO2m_i2(k, n) = sum(states_all{n}{k}.FlowProps.ComponentPhaseMass{2,1}(id_inj2{n}));
    end
end

% Plot
fh = figure(4);
subplot(1,2,1)
hold on
%plot(tsh, disCO2m_i1(:,1)*1e3, '-', 'color', [0.7 0.7 0.7], ...
%     'linewidth', 2, 'DisplayName', '$h=1$ [mm]')
plot(tsh, disCO2m_i1(:,2)*1e3, '-k', 'linewidth', 1.5, 'DisplayName', '$h=2$')
plot(tsh, disCO2m_i1(:,3)*1e3, '-r', 'linewidth', 1, 'DisplayName', '$h=4$')
plot(tsh, disCO2m_i1(:,4)*1e3, '-b', 'linewidth', 0.5, 'DisplayName', '$h=8$')
xlabel('$t$ [h]', latx{:}, 'fontsize', 14)
ylabel('dissolved CO$_2$ [g]', latx{:}, 'fontSize', 14)
%xlim([0 10]), xticks(0:10)
%set(gca,'xscale','log')
legend(latx{:}, 'fontsize', 11)
grid on
hold off

subplot(1,2,2)
hold on
%plot(tsh, disCO2m_i2(:,1)*1e3, '-', 'color', [0.7 0.7 0.7], ...
%     'linewidth', 2, 'DisplayName', '$h=1$ [mm]')
plot(tsh, disCO2m_i2(:,2)*1e3, '-k', 'linewidth', 1.5, 'DisplayName', '$h=2$')
plot(tsh, disCO2m_i2(:,3)*1e3, '-r', 'linewidth', 1, 'DisplayName', '$h=4$')
plot(tsh, disCO2m_i2(:,4)*1e3, '-b', 'linewidth', 0.5, 'DisplayName', '$h=8$')
xlabel('$t$ [h]', latx{:}, 'fontsize', 14)
ylabel('dissolved CO$_2$ [g]', latx{:}, 'fontSize', 14)
%xlim([0 10]), xticks(0:20)
%set(gca,'xscale','log')
legend(latx{:}, 'fontsize', 11)
grid on
hold off
set(fh, 'position', [200, 200, 600, 300]);
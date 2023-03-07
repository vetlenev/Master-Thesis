%% Time plots of states
%
%
clc, clear, close all force 
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers
       
%% Read and organize data
% 1. load data for simulation x
dirs = {'F:\Lluis\mrstSimData\gomBasinFaultLeak\sc1\';
        'F:\Lluis\mrstSimData\gomBasinFaultLeak\sc2\baseCase_nextToF\';
        'F:\Lluis\mrstSimData\gomBasinFaultLeak\sc2\predict\agu22\';
        'F:\Lluis\mrstSimData\gomBasinFaultLeak\sc2\predict\agu22\'};
%fname = {'sc1_NextToF_sc2Mesh63x23523_SGRauto_Spe02_krBase_pvtBase_30y200y'; ...
%         'sc2_baseCase_nextToF_sc2Mesh63x23523_SGRauto_Spe02_krHyst_pvtBase_30y200y'};
problemName = {'proMod_sc1_NextToF_sc2Mesh63x23523_SGRauto_Spe02_80C_krHyst_50y500y'; ...
               'proMod_sc2_baseCase_nextToF_sc2Mesh63x23523_SGRauto_Spe02_krHyst_pvtBase_30y200y'; ...
               'proMod_sc2_predict_nextToF_sc2Mesh63x23523_highPerm_noHyst_pvtBase_30y200y';
               'proMod_sc2_predict_nextToF_sc2Mesh63x23523_lowPerm_noHyst_pvtBase_30y200y'};
inExternalDrive = [1 1 1 1];
numSim = numel(problemName);
resVar = cell(numSim, 1);
for n = 1:numSim
    data = load([dirs{n} problemName{n}]);
    if inExternalDrive(n)
        data.problem.OutputHandlers.states.dataDirectory = dirs{n};
        data.problem.OutputHandlers.reports.dataDirectory = dirs{n};
        data.problem.OutputHandlers.wellSols.dataDirectory = dirs{n};
    end
    [~, states, ~] = getPackedSimulatorOutput(data.problem);
    sg = zeros(data.model.G.cells.num, numel(data.timesteps));
    rs = sg; dp = sg; gmass = sg; pcap = sg;
    for j = 1:numel(data.timesteps)
        sg(:, j) = states{j}.s(:, 2);
        rs(:, j) = states{j}.rs;
        dp(:, j) = states{j}.pressure - data.state0.pressure;
        gmass(:, j) = states{j}.FlowProps.ComponentTotalMass{2};
        pc = data.model.getProp(states{j}, 'CapillaryPressure'); 
        pcap(:, j) = pc{end}; 
    end
    resVar{n}.tstep = cumsum(data.timesteps)/year;
    %resVar{n}.t = t;
    resVar{n}.fluid = data.model.fluid;
    resVar{n}.sg = sg;
    resVar{n}.rs = rs;
    resVar{n}.dp = dp;
    resVar{n}.gmass = gmass;
    resVar{n}.pc = pcap;
    resVar{n}.G = data.model.G;
    %resVar{n}.fault = fault;
    resVar{n}.state0 = data.state0;
    %resVar{nSim}.studyCase = studyCase;
    resVar{n}.W = data.W;
    resVar{n}.ucids = data.ucids;
    
    co2massTot = sum(gmass(:,end));
    disp(['Total Inj CO2 mass = ' num2str(co2massTot/10^9) ' Mt']);

% DON'T FORGET TO CLEAR EVERYTHING BUT THE VARIABLES TO PLOT BEFORE LOADING NEXT SIM, EACH ONE TAKES
% SEVERAL GB OF RAM
clear states sg rs dp gmass pc pcap fault figs data 
disp(['Loading data for simulation ' num2str(n) '/' num2str(numSim) ' finished.'])
end

ucids = resVar{2}.ucids;

%% Plots
xin = [0, 50, 50, 0 , 0]; yin = [0, 0, 6, 6, 0];
colors = [150 150 150; 49 147 206; 92 64 146; 163 31 52]/255;
lst = '-';
clrs = turbo(numSim);
mst = {'d', 's'};
msz = [8, 8];
latx  = {'Interpreter','latex'};
tid = [44, 56];
zcont = [1974, 1871, 1450, 1355, 547.9]/1000;
pctMassCO2 = cell(numSim, 1);
for n = 1:numSim
   % Compute
   % [resRefLeft, resOverlying1, resOverlying2, MMUM, fault]
   massCO2 = [sum(resVar{n}.gmass([ucids.unit_cell_ids{34}, ucids.unit_cell_ids{36}], :))/10^9; ...
              sum(resVar{n}.gmass(ucids.unit_cell_ids{3}, :))/10^9; ...
              sum(resVar{n}.gmass(ucids.unit_cell_ids{5}, :))/10^9; ...
              sum(resVar{n}.gmass(ucids.unit_cell_ids{7}, :))/10^9; ...
              sum(resVar{n}.gmass(ucids.unit_cell_ids{23}, :))/10^9; ...
              sum(resVar{n}.gmass(ucids.unit_cell_ids{61}, :))/10^9]; % [Mt]              
   pctMassCO2{n} = (massCO2 ./ 30)*100;   % 30 Mt of injected CO2.
   
   % Plot
   if n==1
        h = figure(1);
        tiledlayout(2, 3, 'Padding', 'tight', 'TileSpacing', 'tight');
   end
   nexttile(1) % -------------
   hold on
   if n==1
    fill(xin, yin, [193 221 245]/255, 'EdgeColor', 'none')
   end
   plot(resVar{n}.tstep, pctMassCO2{n}(1, :), 'lineStyle', lst, 'color', clrs(n,:), 'lineWidth', 1.5)
   grid on
   ax = gca; ax.YAxis.FontSize = 10; ax.XAxis.FontSize = 10;
   xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
   ylabel('$\% \, M_\mathrm{inj}^{\mathrm{CO}_2}$','fontsize', 12, latx{:})
   %ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
   set(ax,'Xtick', 0:100:500); xlim([0 200])
   set(ax,'Ytick', 0:6)
   ylim([0, max(yin)])
   title('SR (FW)', 'fontweight', 'normal', 'fontsize', 12);
   if n == numSim
       legend('$t_\mathrm{inj}$', 'S$_1$', 'S$_2$ (base)', ...
              'S$_2$ (P$_{\mathrm{high}}$)', 'S$_2$ (P$_{\mathrm{low}}$)', ...
              'interpreter', 'latex', 'location', 'best', 'fontsize', 10)
%        legend('t$_\mathrm{inj}$', 'SR (FW)', 'OR$_1$', 'OR$_2$', 'MMUM', 'Fault', latx{:}, ...
%               'location', 'east', 'box', 'off', 'fontsize', 12, ...
%               'textcolor', [100 100 100]/255)
%       text(20, 4, 'case 2', 'FontSize', 12, 'color', [100 100 100]/255, latx{:})
       hold off
   end
   
   nexttile(2) % -------------
   hold on
   if n==1
    fill(xin, yin, [193 221 245]/255, 'EdgeColor', 'none')
   end
   plot(resVar{n}.tstep, pctMassCO2{n}(2, :), 'lineStyle', lst, 'color', clrs(n,:), 'lineWidth', 1.5)
   grid on
   ax = gca; ax.YAxis.FontSize = 10; ax.XAxis.FontSize = 10;
   %xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
   %ylabel('$\% \, M_\mathrm{inj}^{\mathrm{CO}_2}$','fontsize', 12, latx{:})
   %ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
   set(ax,'Xtick', 0:50:500); xlim([0 200])
   set(ax,'Ytick', 0:0.1:0.5)
   ylim([0, 0.5])
   title('OR_1', 'fontweight', 'normal', 'fontsize', 12);
   
   nexttile(3) % -------------
   hold on
   if n==1
    fill(xin, yin, [193 221 245]/255, 'EdgeColor', 'none')
   end
   plot(resVar{n}.tstep, pctMassCO2{n}(3, :), 'lineStyle', lst, 'color', clrs(n,:), 'lineWidth', 1.5)
   grid on
   ax = gca; ax.YAxis.FontSize = 10; ax.XAxis.FontSize = 10;
   %xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
   %ylabel('$\% \, M_\mathrm{inj}^{\mathrm{CO}_2}$','fontsize', 12, latx{:})
   %ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
   set(ax,'Xtick', 0:50:500); xlim([0 200])
   set(ax,'Ytick', 0:0.1:0.5)
   ylim([0, 0.5])
   title('OR_2', 'fontweight', 'normal', 'fontsize', 12);
   
   nexttile(4) % -------------
   hold on
   if n==1
    fill(xin, yin, [193 221 245]/255, 'EdgeColor', 'none')
   end
   plot(resVar{n}.tstep, pctMassCO2{n}(4, :), 'lineStyle', lst, 'color', clrs(n,:), 'lineWidth', 1.5)
   grid on
   ax = gca; ax.YAxis.FontSize = 10; ax.XAxis.FontSize = 10;
   %xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
   %ylabel('$\% \, M_\mathrm{inj}^{\mathrm{CO}_2}$','fontsize', 12, latx{:})
   %ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
   set(ax,'Xtick', 0:50:500); xlim([0 200])
   set(ax,'Ytick', 0:0.1:0.5)
   ylim([0, 0.5])
   title('OR_3', 'fontweight', 'normal', 'fontsize', 12);
   
   nexttile(5) % -------------
   hold on
   if n==1
    fill(xin, yin, [193 221 245]/255, 'EdgeColor', 'none')
   end
   plot(resVar{n}.tstep, pctMassCO2{n}(5, :), 'lineStyle', lst, 'color', clrs(n,:), 'lineWidth', 1.5)
   grid on
   ax = gca; ax.YAxis.FontSize = 10; ax.XAxis.FontSize = 10;
   %xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
   %ylabel('$\% \, M_\mathrm{inj}^{\mathrm{CO}_2}$','fontsize', 12, latx{:})
   %ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
   set(ax,'Xtick', 0:50:500); xlim([0 200])
   set(ax,'Ytick', 0:0.1:0.5)
   ylim([0, 0.5])
   title('MMUM', 'fontweight', 'normal', 'fontsize', 12);
   
   nexttile(6) % -------------
   hold on
   if n==1
    fill(xin, yin, [193 221 245]/255, 'EdgeColor', 'none')
   end
   plot(resVar{n}.tstep, pctMassCO2{n}(6, :), 'lineStyle', lst, 'color', clrs(n,:), 'lineWidth', 1.5)
   grid on
   ax = gca; ax.YAxis.FontSize = 10; ax.XAxis.FontSize = 10;
   %xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
   %ylabel('$\% \, M_\mathrm{inj}^{\mathrm{CO}_2}$','fontsize', 12, latx{:})
   %ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
   set(ax,'Xtick', 0:50:500); xlim([0 200])
   set(ax,'Ytick', 0:0.1:0.5)
   ylim([0, 0.5])
   title('Fault', 'fontweight', 'normal', 'fontsize', 12);
   set(h, 'Position', [50, 50, 600, 400])
   
   
%    h2 = figure(2);
%    if n == 1
%       p1 = plot([0 0.7], [zcont(1) zcont(1)], '-', 'color', colors(1, :), 'LineWidth', 1.5);
%       hold on 
%       p2 = plot([0 0.7], [zcont(2) zcont(2)], '-k', 'LineWidth', 1.5);
%       plot([0 0.7], [zcont(3) zcont(3)], '-', 'color', colors(1, :), 'LineWidth', 1.5)
%       plot([0 0.7], [zcont(4) zcont(4)], '-k', 'LineWidth', 1.5)
%       plot([0 0.7], [zcont(5) zcont(5)], '-', 'color', colors(1, :), 'LineWidth', 1.5)
%    end
%    fc = resVar{n}.fault.fcells;
%    p3(n) = plot(resVar{n}.sg(fc,tid(n)), resVar{n}.G.cells.centroids(fc, 3)/1000, 'Marker', mst{n}, ...
%         'MarkerSize', msz(n), 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none');
%    
%    p4(n) = plot(resVar{n}.sg(fc,end), resVar{n}.G.cells.centroids(fc, 3)/1000, 'Marker', mst{n}, ...
%         'MarkerSize', msz(n), 'MarkerFaceColor', colors(3,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none');
%    grid on
%    set(gca, 'Ydir', 'reverse')
%    ax = gca; ax.YAxis.FontSize = 32; ax.XAxis.FontSize = 32;
%    xlabel('S$_\mathrm{g}$ [-]', latx{:}, 'fontsize', 38)
%    ylabel('F Cell Depth [km]','fontsize', 38, latx{:})
%    ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
%    set(ax,'Xtick',[0.1 0.3 0.5 0.7])
%    set(ax,'Ytick',[1 1.5 2 2.5])
%    ylim([1 2.5])
%    xlim([ 0 0.7])
%    if n == 2
%        legend([p3(1) p4(1) p3(2) ], 'case 1: t$_\mathrm{inj}$', 'case 1: t$_\mathrm{sim}$', 'case 2', latx{:}, ...
%               'location', 'north', 'box', 'on', 'color', 'w', 'EdgeColor', 'w', 'fontsize', 28, ...
%               'textcolor', [100 100 100]/255)
%    end
%    set(h2, 'Position', [600, 600, 500, 500])
end
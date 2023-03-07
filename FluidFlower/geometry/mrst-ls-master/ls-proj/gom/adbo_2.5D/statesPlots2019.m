%% Time plots of states
%
%
clc, clear, close all force 
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers
       
%% Read and organize data
% 1. load data for simulation x
fid = {'C:\Users\lsalo\matlab\sim_data\scenario_1\mrst2019a_gom_adbo_2.5DTriMesh_SC1_CASE1_1mty_spe02perm_wideX_pvMult_20yInjTime_100ySimTime\spe02_38x45_case1.mat';
       'C:\Users\lsalo\matlab\sim_data\scenario_1\mrst2019a_gom_adbo_2.5DTriMesh_SC1_CASE2_1mty_pcFaultScaled_spe02perm_wideX_pvMult_20yInjTime_100ySimTime\case2_38x45.mat'};
numSim = numel(fid);
resVar = cell(numSim, 1);
for nSim = 1:numSim
    load(fid{nSim})
    for j = 1:numel(timesteps)
        sg(:, j) = states{j}.s(:, 2);
        rs(:, j) = states{j}.rs;
        dp(:, j) = states{j}.dp;
        gmass(:, j) = states{j}.FlowProps.ComponentTotalMass{2};
        pcap(:, j) = states{j}.pc;        
    end
    resVar{nSim}.tstep = cumsum(timesteps)/year;
    resVar{nSim}.t = t;
    resVar{nSim}.fluid = fluid;
    resVar{nSim}.sg = sg;
    resVar{nSim}.rs = rs;
    resVar{nSim}.dp = dp;
    resVar{nSim}.gmass = gmass;
    resVar{nSim}.pc = pcap;
    resVar{nSim}.G = G;
    resVar{nSim}.fault = fault;
    resVar{nSim}.ucids = ucids;
    resVar{nSim}.state0 = state0;
    %resVar{nSim}.studyCase = studyCase;
    resVar{nSim}.W = W;
    
    if nSim == 1
        resLM2Left = load('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/extr_tri/ucids_resLM2Left.mat');
        resVar{nSim}.LM2Left = resLM2Left.resLM2Left;
        resVar{nSim}.mmumId = 6;
    elseif nSim == 2
        resLM2Left = load('mrst-2019a/myprojects/gom/adbo_2.5D/2.5Dmesh/extr_tri/ucids_resLM2LeftUpper.mat');
        resVar{nSim}.LM2Left = resLM2Left.resLM2Left;
        resVar{nSim}.mmumId = 3;
    end

% DON'T FORGET TO CLEAR EVERYTHING BUT THE VARIABLES TO PLOT BEFORE LOADING NEXT SIM, EACH ONE TAKES
% SEVERAL GB OF RAM
clear states sg rs dp gmass pcap fault figs fluid fname G mesh model nls report rock schedule wells wellSols
end

%% Plots
xin = [0, 20, 20, 0 , 0]; yin = [0, 0, 5, 5, 0];
colors = [150 150 150; 49 147 206; 92 64 146; 163 31 52]/255;
lst = {'-', '--'};
mst = {'d', 's'};
msz = [8, 8];
latx  = {'Interpreter','latex'};
tid = [44, 56];
zcont = [1974, 1871, 1450, 1355, 547.9]/1000;
for n = 1:nSim
   % Compute
   CO2MassLeft = sum(resVar{n}.gmass(resVar{n}.LM2Left, :))/10^9;
   CO2MassFault = sum(resVar{n}.gmass(resVar{n}.fault.fcells, :))/10^9;
   CO2MassMMUM = sum(resVar{n}.gmass(resVar{n}.ucids.unit_cell_ids{resVar{n}.mmumId}, :))/10^9;
   %CO2Mass = sum(resVar{n}.gmass)/10^6;
   
   % Plot
   h = figure(1);
   if n == 1
        fill(xin, yin, [193 221 245]/255, 'EdgeColor', 'none')
        hold on
   end
   %plot(resVar{n}.tstep, CO2Mass, 'lineStyle', lst{n}, 'color', colors(1,:), 'lineWidth', 1.5)
   plot(resVar{n}.tstep, CO2MassLeft, 'lineStyle', lst{n}, 'color', colors(2,:), 'lineWidth', 1.5)
   plot(resVar{n}.tstep, CO2MassFault, 'lineStyle', lst{n}, 'color', colors(3,:), 'lineWidth', 1.5)
   plot(resVar{n}.tstep, CO2MassMMUM, 'lineStyle', lst{n}, 'color', colors(4,:), 'lineWidth', 1.5)
   grid on
   ax = gca; ax.YAxis.FontSize = 32; ax.XAxis.FontSize = 32;
   xlabel('t [y]', latx{:}, 'fontsize', 38)
   ylabel('CO$_2$ mass [Mt]','fontsize', 38, latx{:})
   ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
   set(ax,'Xtick',[0 20 40 60 80 100])
   set(ax,'Ytick',[0 1 2 3 4 5])
   ylim([0, max(yin)])
   if n == 2
       legend('t$_\mathrm{inj}$', 'SR (FW)', 'Fault', 'Upper', latx{:}, ...
              'location', 'north', 'box', 'off', 'fontsize', 28, ...
              'textcolor', [100 100 100]/255)
       text(20, 4, 'case 2', 'FontSize', 28, 'color', [100 100 100]/255, latx{:})
   end
   set(h, 'Position', [600, 600, 500, 500])
   
   
   h2 = figure(2);
   if n == 1
      p1 = plot([0 0.7], [zcont(1) zcont(1)], '-', 'color', colors(1, :), 'LineWidth', 1.5);
      hold on 
      p2 = plot([0 0.7], [zcont(2) zcont(2)], '-k', 'LineWidth', 1.5);
      plot([0 0.7], [zcont(3) zcont(3)], '-', 'color', colors(1, :), 'LineWidth', 1.5)
      plot([0 0.7], [zcont(4) zcont(4)], '-k', 'LineWidth', 1.5)
      plot([0 0.7], [zcont(5) zcont(5)], '-', 'color', colors(1, :), 'LineWidth', 1.5)
   end
   fc = resVar{n}.fault.fcells;
   p3(n) = plot(resVar{n}.sg(fc,tid(n)), resVar{n}.G.cells.centroids(fc, 3)/1000, 'Marker', mst{n}, ...
        'MarkerSize', msz(n), 'MarkerFaceColor', colors(2,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none');
   
   p4(n) = plot(resVar{n}.sg(fc,end), resVar{n}.G.cells.centroids(fc, 3)/1000, 'Marker', mst{n}, ...
        'MarkerSize', msz(n), 'MarkerFaceColor', colors(3,:), 'MarkerEdgeColor', 'none', 'LineStyle', 'none');
   grid on
   set(gca, 'Ydir', 'reverse')
   ax = gca; ax.YAxis.FontSize = 32; ax.XAxis.FontSize = 32;
   xlabel('S$_\mathrm{g}$ [-]', latx{:}, 'fontsize', 38)
   ylabel('F Cell Depth [km]','fontsize', 38, latx{:})
   ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
   set(ax,'Xtick',[0.1 0.3 0.5 0.7])
   set(ax,'Ytick',[1 1.5 2 2.5])
   ylim([1 2.5])
   xlim([ 0 0.7])
   if n == 2
       legend([p3(1) p4(1) p3(2) ], 'case 1: t$_\mathrm{inj}$', 'case 1: t$_\mathrm{sim}$', 'case 2', latx{:}, ...
              'location', 'north', 'box', 'on', 'color', 'w', 'EdgeColor', 'w', 'fontsize', 28, ...
              'textcolor', [100 100 100]/255)
   end
   set(h2, 'Position', [600, 600, 500, 500])
end
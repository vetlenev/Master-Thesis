function monitoringCells = plotLine_bo(G, timesteps, W, states, state0, ...
                                       bhp, watchedLoc, ucids, studyCase)
%
%

latx = {'Interpreter', 'latex'};
reportTimes = cumsum(timesteps);

%% BHP
figure(15)
subplot(2,1,1)
plot([0 reportTimes(1:end)/year], [state0.pressure(W.cells)/10^5 bhp(1:end)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 3)
hold on
text(20, 250, 'EoI', latx{:}, 'fontSize', 11)
hold off
xlabel('$t$ [y]', latx{:})
ylabel('$p$ [bar]', latx{:})
%ylim([100 171])
%yticks([100 110 120 130 140 150 160 170])
grid on
ylim([200 400])
title('BHP vs time', latx{:}, 'fontSize', 12)

subplot(2,1,2)
plot([0 reportTimes(1:15)/day], [state0.pressure(W.cells)/10^5 bhp(1:15)./10^5], '-ok', 'markerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 3)
xlabel('$t$ [day]', latx{:})
xticks([0 7 30 60 90 120 150 180 210 240 270 300 330 365])
%yticks([100 110 120 130 140 150 160 165])
%ylim([100 165])
hold on
text(120, 380, 'Detail during first year', latx{:}, 'fontSize', 11)
hold off
ylabel('$p$ [bar]', latx{:})
grid on
xlim([0 365])
ylim([200 400])

%% Dp vs time at selected cells
dist = pdist2(G.cells.centroids, watchedLoc);
[~, monitoringCells] = min(dist,[],1);      
monitoringCells = [W.cells monitoringCells];

dpCells = zeros(numel(states), numel(monitoringCells));
sCells = zeros(numel(states), numel(monitoringCells));
c_refr =  ucids.unit_cell_ids{35}(~ismember(ucids.unit_cell_ids{35}, ...
                                            ucids.unit_cell_ids{36})); %grid sc2!
dp_res = zeros(numel(states), 1);
sg_res = zeros(numel(states), 1);
for n=1:numel(states)
    dpCells(n, :) = states{n}.dp(monitoringCells);
    sCells(n, :) = states{n}.s(monitoringCells, 2);
    if strcmp(studyCase, 'nextToF_FW')
        dp_res(n) = mean(states{n}.dp(ucids.unit_cell_ids{36})); % grid sc2!
        sg_res(n) = mean(states{n}.s(ucids.unit_cell_ids{36}, 2));
    elseif strcmp(studyCase, 'nextToF')
        dp_res(n) = mean(states{n}.dp(c_refr));
        sg_res(n) = mean(states{n}.s(c_refr, 2));
    end
end

f16 = figure(6);
left_color = [0 0 0];
right_color = [0 0 1];
set(f16,'defaultAxesColorOrder',[left_color; right_color]);
tiledlayout(1, 5, 'Padding', 'tight', 'TileSpacing', 'tight');
nexttile
fill([-5 -5 50 50 -5], [-1 1.2 1.2 -1 -1], [0 0.4470 0.7410], 'faceAlpha', 0.25, 'EdgeColor', 'none')
hold on
yyaxis left
plot([0 reportTimes(1:end)/year], [0; dp_res./10^5], '-ok', 'markerFaceColor', [0.7 0.7 0.7], 'markerSize', 5)
%annotation(f16,'doublearrow', [0.130697210391006 0.157422969187675], [0.282333333333333 0.280952380952381]);
%text(5, 0, 'Inj.', latx{:}, 'fontSize', 12)
hold off
xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
xticks([0 50 115 200:100:300])
xlim([-5 300])
ylim([0 1.2])
ylabel('$\Delta p$ [bar]', latx{:}, 'fontsize', 12)
yyaxis right
plot([0 reportTimes(1:end)/year], [0; sg_res], '-sb', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 4)
ylabel('$S_\mathrm{g}$ [-]', latx{:}, 'fontsize', 12)
ylim([0 0.7])
grid on
title('Avg SR', 'fontweight', 'normal', 'fontSize', 11)

nexttile
fill([-5 -5 50 50 -5], [-5 30 30 -5 -5], [0 0.4470 0.7410], 'faceAlpha', 0.25, 'EdgeColor', 'none')
hold on
yyaxis left
plot([0 reportTimes(1:end)/year], [0; dpCells(1:end,1)./10^5], '-ok', 'markerFaceColor', [0.7 0.7 0.7], 'markerSize', 5)
%annotation(f16,'doublearrow', [0.130697210391006 0.157422969187675], [0.282333333333333 0.280952380952381]);
%text(5, 0, 'Inj.', latx{:}, 'fontSize', 12)
hold off
%xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
xticks([0 50 115 200:100:500])
xlim([-5 300])
ylim([-5 30])
%ylabel('$\Delta p$ [bar]', latx{:}, 'fontsize', 12)
yyaxis right
plot([0 reportTimes(1:end)/year], [0; sCells(:, 1)], '-sb', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 4)
%ylabel('$S_\mathrm{g}$ [-]', latx{:}, 'fontsize', 12)
ylim([0 0.7])
grid on
title('Injector', 'fontweight', 'normal', 'fontSize', 11)

nexttile
fill([-5 -5 50 50 -5], [-5 25 25 -5 -5], [0 0.4470 0.7410], 'faceAlpha', 0.25, 'EdgeColor', 'none')
hold on
yyaxis left
plot([0 reportTimes(1:end)/year], [0; dpCells(1:end,2)./10^5], '-ok', 'markerFaceColor', [0.7 0.7 0.7], 'markerSize', 5)
%xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
%ylabel('$\Delta p$ [bar]', latx{:}, 'fontsize', 12)
ylim([0 8])
xticks([0 50 115 200:100:500])
xlim([-5 300])
yyaxis right
plot([0 reportTimes(1:end)/year], [0; sCells(:, 2)], '-sb', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 4)
%ylabel('$S_\mathrm{g}$ [-]', latx{:}, 'fontsize', 12)
ylim([0 0.7])
grid on
hold off
%legend('9.7m', '97m', '1013m', '5014m', '24.834km', '49.684km')
%ylim([0 100])
title('Top SR', 'fontweight', 'normal', 'fontSize', 11)

nexttile
fill([-5 -5 50 50 -5], [-5 25 25 -5 -5], [0 0.4470 0.7410], 'faceAlpha', 0.25, 'EdgeColor', 'none')
hold on
yyaxis left
plot([0 reportTimes(1:end)/year], [0; dpCells(1:end,3)./10^5], '-ok', 'markerFaceColor', [0.7 0.7 0.7], 'markerSize', 5)
%xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
%ylabel('$\Delta p$ [bar]', latx{:}, 'fontsize', 12)
ylim([0 8])
xticks([0 50 115 200:100:500])
xlim([-5 300])
yyaxis right
plot([0 reportTimes(1:end)/year], [0; sCells(:, 3)], '-sb', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 4)
%ylabel('$S_\mathrm{g}$ [-]', latx{:}, 'fontsize', 12)
ylim([0 0.7])
grid on
hold off
%legend('9.7m', '97m', '1013m', '5014m', '24.834km', '49.684km')
%ylim([0 100])
title('Base TS', 'fontweight', 'normal', 'fontSize', 11)

nexttile
fill([-5 -5 50 50 -5], [-5 25 25 -5 -5], [0 0.4470 0.7410], 'faceAlpha', 0.25, 'EdgeColor', 'none')
hold on
yyaxis left
plot([0 reportTimes(1:end)/year], [0; dpCells(1:end,4)./10^5], '-ok', 'markerFaceColor', [0.7 0.7 0.7], 'markerSize', 5)
%xlabel('$t$ [y]', latx{:}, 'fontsize', 12)
%ylabel('$\Delta p$ [bar]', latx{:}, 'fontsize', 12)
ylim([-1 8])
xticks([0 50 115 200:100:500])
xlim([-5 300])
yyaxis right
plot([0 reportTimes(1:end)/year], [0; sCells(:, 4)], '-sb', 'markerFaceColor', [125, 255, 255]/255, 'markerSize', 4)
%ylabel('$S_\mathrm{g}$ [-]', latx{:}, 'fontsize', 12)
ylim([0 0.7])
grid on
hold off
%legend('9.7m', '97m', '1013m', '5014m', '24.834km', '49.684km')
%ylim([0 100])
title('Fault (SR)', 'fontweight', 'normal', 'fontSize', 11)
set(f16,'position',[50 50 1000 250],...
    'paperunits','points','papersize', [1000 250])

end
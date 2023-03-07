%
% Plot mass
% Fig. 1: plot mass in whole domain
% Fig. 2: mass in box A and box b
%
clear, close all

% Simulation data
pth = 'C:/Users/lsalo/matlab/sim_data/mrst/fluidflower/results/benchmark/';
dirnames = {'bmesh5_thickVar_removeS1_modelcase1_D1e-09_Im1_SgMinHyst0.02_kdF3mm_pcD0.3/';
            'bmesh5_thickVar_removeS1_modelcase2_D1e-09_Im1_SgMinHyst0.02/';
            'bmesh5_thickVar_removeS1_modelcase3_D1e-09_Im1_SgMinHyst0.02/';
            'bmesh5_thickVar_removeS1_modelcase3_D3e-09_Im1_SgMinHyst0.02/'};
fname = {'mass', 'time_series'};

for n=1:2
    if n==1 % total mass
        mass1 = readmatrix([pth dirnames{1} fname{1} '.csv']);
        mass2 = readmatrix([pth dirnames{2} fname{1} '.csv']);
        mass3 = readmatrix([pth dirnames{3} fname{1} '.csv']);
        mass33 = readmatrix([pth dirnames{4} fname{1} '.csv']);
    elseif n==2 % Box A and Box B
        tseries1 = readmatrix([pth dirnames{1} fname{2} '.csv']);
        tseries2 = readmatrix([pth dirnames{2} fname{2} '.csv']);
        tseries3 = readmatrix([pth dirnames{3} fname{2} '.csv']);
        tseries33 = readmatrix([pth dirnames{4} fname{2} '.csv']);
    end
end

% Plots
% --------------------------- Total Mass ---------------------------------
latx = {'interpreter', 'latex'};
h = figure(35);
subplot(1,4,1)
hold on
p1 = plot(mass1(:,1)/hour, mass1(:,2)*1e3, '-b', ...
          'displayname', 'm$_\mathrm{I}$');
p2 = plot(mass2(:,1)/hour, mass2(:,2)*1e3, '-c', ...
          'displayname', 'm$_\mathrm{II}$'); 
p3 = plot(mass3(:,1)/hour, mass3(:,2)*1e3, '-m', ...
          'displayname', 'm$_{\mathrm{III},1}$');
p4 = plot(mass33(:,1)/hour, mass33(:,2)*1e3, '--m', ...
          'displayname', 'm$_{\mathrm{III},3}$');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(12, 0.2, 'Mobile CO$_2$', latx{:}, 'fontsize', 11)
legend([p1 p2 p3 p4], latx{:}, 'fontsize', 10)

subplot(1,4,2)
hold on
plot(mass1(:,1)/hour, mass1(:,3)*1e3, '-b');
plot(mass2(:,1)/hour, mass2(:,3)*1e3, '-c'); 
plot(mass3(:,1)/hour, mass3(:,3)*1e3, '-m');
plot(mass33(:,1)/hour, mass33(:,3)*1e3, '--m');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
ylim([0 0.1]); yticks(0:.01:0.1)
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(24, 0.09, 'Immobile CO$_2$', latx{:}, 'fontsize', 11)

subplot(1,4,3)
hold on
plot(mass1(:,1)/hour, mass1(:,4)*1e3, '-b');
plot(mass2(:,1)/hour, mass2(:,4)*1e3, '-c'); 
plot(mass3(:,1)/hour, mass3(:,4)*1e3, '-m');
plot(mass33(:,1)/hour, mass33(:,4)*1e3, '--m');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
ylim([0 9]); yticks([0:2:8 9])
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(12, 0.5, 'Dissolved CO$_2$', latx{:}, 'fontsize', 11)

subplot(1,4,4)
hold on
plot(mass1(:,1)/hour, mass1(:,5)*1e3, '-b');
plot(mass2(:,1)/hour, mass2(:,5)*1e3, '-c'); 
plot(mass3(:,1)/hour, mass3(:,5)*1e3, '-m');
plot(mass33(:,1)/hour, mass33(:,5)*1e3, '--m');
grid on
xlim([0 120]); xticks([0 24 48 72 96 120]);
ylim([0 9]); yticks([0:2:7 8.13 9])
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(24, 0.5, 'Total CO$_2$', latx{:}, 'fontsize', 11)
set(h,'position',[100 100 1200 250])


% -------------------------- Mass in A and B -----------------------------
totCO2A1 = tseries1(:,4) + tseries1(:,5) + tseries1(:,6);
totCO2A2 = tseries2(:,4) + tseries2(:,5) + tseries2(:,6);
totCO2A3 = tseries3(:,4) + tseries3(:,5) + tseries3(:,6);
totCO2A33 = tseries33(:,4) + tseries33(:,5) + tseries33(:,6);
totCO2B1 = tseries1(:,8) + tseries1(:,9) + tseries1(:,10);
totCO2B2 = tseries2(:,8) + tseries2(:,9) + tseries2(:,10);
totCO2B3 = tseries3(:,8) + tseries3(:,9) + tseries3(:,10);
totCO2B33 = tseries33(:,8) + tseries33(:,9) + tseries33(:,10);
h = figure(36);
subplot(2,4,1)
hold on
p1 = plot(tseries1(:,1)/hour, tseries1(:,4)*1e3, '-b', ...
          'displayname', 'm$_\mathrm{I}$');
p2 = plot(tseries2(:,1)/hour, tseries2(:,4)*1e3, '-c', ...
          'displayname', 'm$_\mathrm{II}$'); 
p3 = plot(tseries3(:,1)/hour, tseries3(:,4)*1e3, '-m', ...
          'displayname', 'm$_{\mathrm{III},1}$');
p4 = plot(tseries3(:,1)/hour, tseries3(:,4)*1e3, '--m', ...
          'displayname', 'm$_{\mathrm{III},3}$');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(12, 0.2, 'Mobile CO$_2$ A', latx{:}, 'fontsize', 11)
legend([p1 p2 p3 p4], latx{:}, 'fontsize', 10)

subplot(2,4,2)
hold on
plot(tseries1(:,1)/hour, tseries1(:,5)*1e3, '-b');
plot(tseries2(:,1)/hour, tseries2(:,5)*1e3, '-c'); 
plot(tseries3(:,1)/hour, tseries3(:,5)*1e3, '-m');
plot(tseries33(:,1)/hour, tseries33(:,5)*1e3, '--m');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(24, 1.75e-3, 'Immobile CO$_2$ A', latx{:}, 'fontsize', 11)

subplot(2,4,3)
hold on
plot(tseries1(:,1)/hour, tseries1(:,6)*1e3, '-b');
plot(tseries2(:,1)/hour, tseries2(:,6)*1e3, '-c'); 
plot(tseries3(:,1)/hour, tseries3(:,6)*1e3, '-m');
plot(tseries33(:,1)/hour, tseries33(:,6)*1e3, '--m');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(12, 0.3, 'Dissolved CO$_2$ A', latx{:}, 'fontsize', 11)

subplot(2,4,4)
hold on
plot(tseries1(:,1)/hour, totCO2A1*1e3, '-b');
plot(tseries2(:,1)/hour, totCO2A2*1e3, '-c'); 
plot(tseries3(:,1)/hour, totCO2A3*1e3, '-m');
plot(tseries33(:,1)/hour, totCO2A33*1e3, '--m');
grid on
xlim([0 120]); xticks([0 24 48 72 96 120]);
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(24, 0.3, 'Total CO$_2$ A', latx{:}, 'fontsize', 11)

subplot(2,4,5)
hold on
plot(tseries1(:,1)/hour, tseries1(:,8)*1e3, '-b');
p2 = plot(tseries2(:,1)/hour, tseries2(:,8)*1e3, '-c'); 
p3 = plot(tseries3(:,1)/hour, tseries3(:,8)*1e3, '-m');
plot(tseries33(:,1)/hour, tseries33(:,8)*1e3, '--m');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(24, 0.017, 'Mobile CO$_2$ B', latx{:}, 'fontsize', 11)

subplot(2,4,6)
hold on
plot(tseries1(:,1)/hour, tseries1(:,9)*1e3, '-b');
plot(tseries2(:,1)/hour, tseries2(:,9)*1e3, '-c'); 
plot(tseries3(:,1)/hour, tseries3(:,9)*1e3, '-m');
plot(tseries33(:,1)/hour, tseries33(:,9)*1e3, '--m');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(24, 5.1e-3, 'Immobile CO$_2$ B', latx{:}, 'fontsize', 11)

subplot(2,4,7)
hold on
plot(tseries1(:,1)/hour, tseries1(:,10)*1e3, '-b');
plot(tseries2(:,1)/hour, tseries2(:,10)*1e3, '-c'); 
plot(tseries3(:,1)/hour, tseries3(:,10)*1e3, '-m');
plot(tseries33(:,1)/hour, tseries33(:,10)*1e3, '--m');
grid on
xlim([0 72]); xticks([0 5 12 24 48 72]);
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(24, 0.85, 'Dissolved CO$_2$ B', latx{:}, 'fontsize', 11)

subplot(2,4,8)
hold on
plot(tseries1(:,1)/hour, totCO2B1*1e3, '-b');
plot(tseries2(:,1)/hour, totCO2B2*1e3, '-c'); 
plot(tseries3(:,1)/hour, totCO2B3*1e3, '-m');
plot(tseries33(:,1)/hour, totCO2B33*1e3, '--m');
grid on
xlim([0 120]); xticks([0 24 48 72 96 120]);
%xlabel('$t$ [h]', latx{:}, 'fontsize', 12);
%ylabel('Mass [g]', latx{:}, 'fontsize', 12);
text(48, 0.85, 'Total CO$_2$ B', latx{:}, 'fontsize', 11)
set(h,'position',[100 100 1100 500])



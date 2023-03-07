%% Clay content and fault permeability (part 1)
%
% Analysis of fault zone clay content vs fault permeability, comparison
% with SGR and SGR-based fault perm models (Stratigraphies from
% stratiControls.m)
%

clear, close all

% Load data from stratiControls.m
load('sim_data/predict/stratiControls/stratiControls_Nsim5000.mat');

% Organize data
fvcl = [cell2mat(cellfun(@(x) x.Vcl, faults_all{1}, 'UniformOutput', false)) ...
        cell2mat(cellfun(@(x) x.Vcl, faults_all{2}, 'UniformOutput', false)) ...
        cell2mat(cellfun(@(x) x.Vcl, faults_all{3}, 'UniformOutput', false))];
K = {log10(cell2mat(cellfun(@(x) x.Perm, faults_all{1}, 'UniformOutput', false))./(milli*darcy)), ...
     log10(cell2mat(cellfun(@(x) x.Perm, faults_all{2}, 'UniformOutput', false))./(milli*darcy)), ...
     log10(cell2mat(cellfun(@(x) x.Perm, faults_all{3}, 'UniformOutput', false))./(milli*darcy))} ;
nedg = 26;
%edges = [linspace(min(fvcl(:,1)), max(fvcl(:,1)), nbins)', ...
%         linspace(min(fvcl(:,2)), max(fvcl(:,2)), nbins)', ...
%         linspace(min(fvcl(:,3)), max(fvcl(:,3)), nbins)'];
edges = linspace(0.2, 0.6, nedg)';
logMinP = [min(min(K{1})), min(min(K{2})), min(min(K{3}))];
logMaxP = [max(max(K{1})), max(max(K{2})), max(max(K{3}))];
edgesk = [linspace(fix(logMinP(1))-1, fix(logMaxP(1))+1, nedg)', ...
          linspace(fix(logMinP(2))-1, fix(logMaxP(2))+1, nedg)', ...
          linspace(fix(logMinP(3))-1, fix(logMaxP(3))+1, nedg)'];
      
% Section SGR and perm
SGR = [getSGR(thickness{1}, vcl{1}), getSGR(thickness{2}, vcl{2}), ...
       getSGR(thickness{3}, vcl{3})];
% Check SGR
% G = makeFaultGrid(faults_all{3}{1}.MatProps.thick, faults_all{3}{1}.Disp, ...
%                   faults_all{3}{1}.Grid.targetCellDim);
% SGR3 = SGR(:,3); 
% hSGR = linspace(0, faults_all{3}{1}.Throw, numel(SGR3))';
% hSGR_inFault = hSGR/sind(faults_all{3}{1}.Dip);
% hG = G.cells.centroids(1:G.cartDims(1):end, 2);
% SGR3 = interp1(hSGR_inFault, SGR3, hG); 
% SGR3 = repelem(SGR3, G.cartDims(1));
% figure(96)                          
% set(gca, 'colormap', hot);
% plotToolbar(G, SGR3, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
% xlim([0 faults_all{3}{1}.MatProps.thick]); ylim([0 faults_all{3}{1}.Disp]); 
% c = colorbar;                         
                          
SGR = [min(SGR(:, 1)), mean(SGR(:, 1)), max(SGR(:, 1)); ...
       min(SGR(:, 2)), mean(SGR(:, 2)), max(SGR(:, 2)); ...
       min(SGR(:, 3)), mean(SGR(:, 3)), max(SGR(:, 3))];
logk_spe02 = log10([getPermSpe02(SGR([1 4 7]), zf(1), zmax{1}{1}(1)); ...
                    getPermSpe02(SGR([2 5 8]), zf(2), zmax{2}{1}(1)); ...
                    getPermSpe02(SGR([3 6 9]), zf(3), zmax{3}{1}(1))]);
                
                
% GoM data (GoM Atlas, 2017, Ch. 3)
% XRD (Mass fraction assumed \approx Vcl since rho_s_c and rho_s_bulk are 
% probably similar)
z = [10578 10609]*0.3048;  
mcl = [0.27 0.38];
n = [0.0315 0.0932];
logk_gom = log10([0.000145 0.00455]);      % log(mD)

% Red Fault (Discuss, not plotted)
 
 
 % NNS data
 % Spe02 is already plotted
 
% Taranaki data (Childs et al., 2007)
vclt = [0.225 0.51];
logk_tar = [-3.699 0];
   
% SGR histograms
labl = "f$_{V_\mathrm{cl}}$ [-]";
latx = {'Interpreter', 'latex'};
sz = [14, 12];
gg = [195, 224, 202; 67, 110, 81]./255;
pp = [175, 170, 189; 92, 76, 135]./255;
bb = [0.6941    0.7961    0.8902; 0.2431    0.3765    0.5020];
fh = figure(1);
tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
clrs = copper(10);

nexttile(1)
hold on
l1 = fill([SGR(1) SGR(7) SGR(7) SGR(1) SGR(1)], [0.001 0.001 0.499 0.499 0.001], ...
     [255, 230, 196]/255, 'edgecolor', 'none', 'FaceAlpha', 0.8, 'DisplayName', 'SGR');
l2 = plot([SGR(4) SGR(4)], [0 0.5], '-', 'color', clrs(1,:), 'DisplayName', ...
     '$\overline{\mathrm{SGR}}$');
l3 = plot([mcl(1) mcl(1)],[0.001 0.499], 'k--', 'lineWidth', 1.25, ...
           'DisplayName', 'Lu17');
plot([mcl(2) mcl(2)],[0.001 0.499], 'k--', 'lineWidth', 1.25);
%plot(mcl(1), 0.25, 'ok', 'markerFaceColor', bb(2, :))
%plot(mcl(2), 0.25, 'ok', 'markerFaceColor', bb(2, :))
histogram(fvcl(:, 1), edges, 'Normalization', 'probability', ...
          'FaceColor', clrs(3, :), 'FaceAlpha', 1)
%xlabel(labl(1), latx{:}, 'fontSize', sz(2))
xlabel('Clay fraction [-]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
%xlim([min([fvcl(:,1); SGR(1)])-0.01, max([fvcl(:,1); SGR(7)])+0.01])
xticks([0.2 0.3 0.4 0.5 0.6]); xlim([0.2 0.6])
%xlim([0.2 0.65]); xticks([0.2 0.3 0.4 0.5 0.6])
ylim([0 0.5]); yticks(0:.1:.5)
grid on
leg = legend([l1 l2 l3], latx{:}, 'fontSize', 9);
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));

nexttile(3)
hold on
fill([SGR(2) SGR(8) SGR(8) SGR(2) SGR(2)], [0.001 0.001 0.499 0.499 0.001], ...
     [255, 230, 196]/255, 'edgecolor', 'none')
plot([SGR(5) SGR(5)], [0 0.5], '-', 'color', clrs(1,:))
histogram(fvcl(:, 2), edges, 'Normalization', 'probability', ...
          'FaceColor', clrs(5, :), 'FaceAlpha', 1)  
%xlabel(labl(1), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
%xlim([min([fvcl(:,2); SGR(2)])-0.01, max([fvcl(:,2); SGR(8)])+0.01])
xticks([0.2 0.3 0.4 0.5 0.6]); xlim([0.2 0.6])
ylim([0 0.5]); yticks(0:.1:.5)
grid on

nexttile(5)
hold on
l1 = plot([vclt(1) vclt(1)], [0.001 0.499], 'k-.', ...
          'lineWidth', 1.25, 'DisplayName', 'Chi07');
plot([vclt(2) vclt(2)], [0.001 0.499], 'k-.', ...
          'lineWidth', 1.25, 'DisplayName', 'Chi07');
fill([SGR(3) SGR(9) SGR(9) SGR(3) SGR(3)], [0.001 0.001 0.499 0.499 0.001], ...
     [255, 230, 196]/255, 'edgecolor', 'none')
plot([SGR(6) SGR(6)], [0 0.5], '-', 'color', clrs(1,:))
histogram(fvcl(:, 3), edges, 'Normalization', 'probability', ...
          'FaceColor', clrs(7, :), 'FaceAlpha', 1)
%xlabel(labl(1), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
%xlim([min([fvcl(:,3); SGR(3)])-0.01, max([fvcl(:,3); SGR(9)])+0.01])
xticks([0.2 0.3 0.4 0.5 0.6]); xlim([0.2 0.6])
ylim([0 0.5]); yticks(0:.1:.5)
grid on
legend(l1, latx{:}, 'fontSize', 9)


% kxx histogram
labls = ["$\log_{10}(k_{xx}$ [mD])", ...
        "$\log_{10}(k_{yy}$ [mD])", ...
        "$\log_{10}(k_{zz}$ [mD])"];
%clrs = [0, 74, 138; 0, 106, 196; 0, 138, 255]/255;

nexttile(2)
hold on
l1 = fill([logk_spe02(7) logk_spe02(1) logk_spe02(1) logk_spe02(7) logk_spe02(7)], ...
     [0.001 0.001 0.499 0.499 0.001], ...
     [0.85 0.85 0.85], 'edgecolor', 'none', 'DisplayName', 'Spe02');
l2 = plot([logk_spe02(4) logk_spe02(4)], [0 0.5], '-k', ...
          'DisplayName', '$\mathrm{Spe02(\overline{\mathrm{SGR}})}$');
l3 = plot([logk_gom(1) logk_gom(1)], [0.001 0.499], ...
           'k--', 'lineWidth', 1.25);
plot([logk_gom(2) logk_gom(2)], [0.001 0.499], ...
     'k--', 'lineWidth', 1.25);
histogram(K{1}(:, 1), edgesk(:,1), 'Normalization', 'probability', ...
    'FaceColor', [0.3 0.3 0.3], 'FaceAlpha', 1)
xlabel(labls(1), latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(min(logMinP(1), logk_spe02(7)))-1 ...
     fix(max(logMaxP(1), logk_spe02(1)))+1])
%xlim([-7 2]); xticks([-6 -4 -2 0 2]);
ylim([0 0.5]); yticks(0:.1:.5)
grid on
leg = legend([l1 l2], latx{:}, 'fontSize', 9);
%set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));

nexttile(4)
hold on
fill([logk_spe02(8) logk_spe02(2) logk_spe02(2) logk_spe02(8) logk_spe02(8)], ...
     [0.001 0.001 0.499 0.499 0.001], ...
     [0.85 0.85 0.85], 'edgecolor', 'none')
plot([logk_spe02(5) logk_spe02(5)], [0 0.5], '-k')
histogram(K{2}(:, 1), edgesk(:,1), 'Normalization', 'probability', ...
    'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
%xlabel(labls(1), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(min(logMinP(2), logk_spe02(8)))-1 ...
      fix(max(logMaxP(2), logk_spe02(2)))+1])
xlim([-7 2]); xticks([-6 -4 -2 0 2]);
ylim([0 0.5]); yticks(0:.1:.5)
grid on

nexttile(6)
hold on
l1 = plot([logk_tar(1) logk_tar(1)], [0.001 0.499], ...
     'k-.', 'lineWidth', 1.25);
plot([logk_tar(2) logk_tar(2)], [0.001 0.499], ...
     'k-.', 'lineWidth', 1.25);
fill([logk_spe02(9) logk_spe02(3) logk_spe02(3) logk_spe02(9) logk_spe02(9)], ...
     [0.001 0.001 0.499 0.499 0.001], ...
     [0.85 0.85 0.85], 'edgecolor', 'none')
plot([logk_spe02(6) logk_spe02(6)], [0 0.5], '-k')
histogram(K{3}(:, 1), edgesk(:,1), 'Normalization', 'probability', ...
    'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 1)
%xlabel(labls(1), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(min(logMinP(3), logk_spe02(9)))-1 ...
      fix(max(logMaxP(3), logk_spe02(3)))+1])
xlim([-7 2]); xticks([-6 -4 -2 0 2]);
ylim([0 0.5]); yticks(0:.1:.5)
grid on
%legend(l1, latx{:}, 'fontSize', 10)

set(fh, 'position', [200, 200, 400, 500]);
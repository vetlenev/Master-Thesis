%% Plot smear-based fault permeability algorithm convergence results
% Save data individually from faultPermSmearThrow_base. Here we load and
% plot perm as a function of cell resolution (length) in the vertical
% direction, since it's the one controlling the square matrix and grid
% dimensions.

clc, clear, close all force 
mrstModule add ls-proj ls-utils coarsegrid upscaling incomp mimetic mpfa ...
           streamlines mrst-gui

% Load data
% 4 beds
% dir = 'C:\Users\lsalo\matlab\sim_data\fault_perm_algorithm\singleParamVariations\convergence_analysis\4beds_rng20\';
% fname = {'resL20_4beds_fT2fdip60_phi1030_SSFc3.6_LsSmear', ...
%          'resL5_4beds_fT2fdip60_phi1030_SSFc3.6_LsSmear', ...
%          'resL2_4beds_fT2fdip60_phi1030_SSFc3.6_LsSmear', ...
%          'resL1_4beds_fT2fdip60_phi1030_SSFc3.6_LsSmear', ...
%          'resL0.5_4beds_fT2fdip60_phi1030_SSFc3.6_LsSmear', ...
%          'resL0.2_4beds_fT2fdip60_phi1030_SSFc3.6_LsSmear', ...
%          'resL0.1_4beds_fT2fdip60_phi1030_SSFc3.6_LsSmear'};

% 10 beds
dir = 'C:\Users\lsalo\matlab\sim_data\fault_perm_algorithm\singleParamVariations\convergence_analysis\10beds_rng5\';
fname = {'resL7.5_10beds_fT2fdip60_phi1030_SSFc6_LsSmear', ...
         'resL5_10beds_fT2fdip60_phi1030_SSFc6_LsSmear', ...
         'resL2_10beds_fT2fdip60_phi1030_SSFc6_LsSmear', ...
         'resL1_10beds_fT2fdip60_phi1030_SSFc6_LsSmear', ...
         'resL0.5_10beds_fT2fdip60_phi1030_SSFc6_LsSmear', ...
         'resL0.2_10beds_fT2fdip60_phi1030_SSFc6_LsSmear', ...
         'resL0.1_10beds_fT2fdip60_phi1030_SSFc6_LsSmear'};

cellDim = zeros(numel(fname), 2);
upscaledPerms = zeros(numel(fname), 3);
for n=1:numel(fname)
   data = load([dir fname{n} '.mat']);
   cellDim(n,:) = data.G.xzFaceDim;
   upscaledPerms(n,:) = data.Perm;   
end

% Plot permeabilities vs 1/h
hL = cellDim(:,2);
A  = cellDim(:,1).*cellDim(:,2);
latx = {'Interpreter','latex'};
k_md = upscaledPerms./(milli*darcy);

f9 = figure(99);
subplot(1,3,1)
plot(1./hL, k_md(:,1), '-ok', 'MarkerSize', 4)
hold on
plot(1/hL(4), k_md(4,1), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
grid on
xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{xx}$ [mD]', latx{:}, 'fontSize', 12)
xlim([0 12])
%ylim([1.1 1.5]*10^-3)
ylim([0 0.2])
xticks([0.001 1 2 5 10])
xticklabels({'10^{-3}', '1' '2' '5' '10'})
subplot(1,3,2)
plot(1./hL, k_md(:,2), '-or', 'MarkerSize', 4)
hold on
plot(1/hL(4), k_md(4,2), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 6)
hold off
grid on
xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{yy}$ [mD]', latx{:}, 'fontSize', 12)
xlim([0 12])
%ylim([8 18])
ylim([40 80])
xticks([0.001 1 2 5 10])
xticklabels({'10^{-3}', '1' '2' '5' '10'})
subplot(1,3,3)
plot(1./hL, k_md(:,3), '-ob')
hold on
plot(1/hL(4), k_md(4,3), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 6)
hold off
grid on
xlabel('$1/h_\mathrm{L}$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{zz}$ [mD]', latx{:}, 'fontSize', 12)
xlim([0 12])
%ylim([0.01 0.018])
ylim([0 50])
xticks([0.001 1 2 5 10])
xticklabels({'10^{-3}', '1' '2' '5' '10'})
set(f9, 'position', [500, 500, 850, 220]);

f10 = figure(100);
subplot(1,3,1)
plot(1./A, k_md(:,1), '-ok', 'MarkerSize', 4)
grid on
hold on
plot(1/A(4), k_md(4,1), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
hold off
xlabel('$1/A$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{xx}$ [mD]', latx{:}, 'fontSize', 12)
xlim([0 6000])
%ylim([1.1 1.5]*10^-3)
ylim([0 0.2])
xticks([10^-4 1000 6000])
xticklabels({'10^{-4}', '10^{3}', '6·10^3'})
subplot(1,3,2)
plot(1./A, k_md(:,2), '-or', 'MarkerSize', 4)
hold on
plot(1/A(4), k_md(4,2), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 6)
hold off
grid on
xlabel('$1/A$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{yy}$ [mD]', latx{:}, 'fontSize', 12)
%xlim([0 6000])
%ylim([8 18])
ylim([40 80])
xticks([10^-4 1000 6000])
xticklabels({'10^{-4}', '10^{3}', '6·10^3'})
subplot(1,3,3)
plot(1./A, k_md(:,3), '-ob', 'MarkerSize', 4)
hold on
plot(1/A(4), k_md(4,3), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 6)
hold off
grid on
xlabel('$1/A$ [m$^{-1}$]', latx{:}, 'fontSize', 12)
ylabel('$\hat{k}_{zz}$ [mD]', latx{:}, 'fontSize', 12)
%xlim([0 6000])
%ylim([0.01 0.018])
ylim([0 50])
xticks([10^-4 1000 6000])
xticklabels({'10^{-4}', '10^{3}', '6·10^3'})
sgtitle('Upscaled permeability vs cell area')
set(f10, 'position', [500, 500, 850, 220]);
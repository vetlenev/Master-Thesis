clear all;
close all;

warning('off','all')
mrstModule add ad-mechanics ad-core ad-props ad-blackoil vemmech deckformat mrst-gui

%% Summary of the options
clear opt
opt.nonlinearTolerance = 1e-6;
opt.splittingTolerance = 1e-6;
opt.verbose            = false;
opt.splittingVerbose   = false;
% write intro text for each case
writeIntroText = @(opt)(fprintf('\n*** Start new simulation\n* fluid model : %s\n* method : %s\n\n', opt.fluid_model, opt.method));

% Check initialization
% opt.grid = 'uniform';           % 'CR' for cappa-rutqvist grid or 'uniform'
% opt.inj_rate = 0.0;             % 0.02 for actual fluid injected
% opt.fluid_model = 'blackoil-wg'; %% 'blackoil-wg', 'blackoil-ccs', 'oil water', 'water'
% opt.method = 'fully coupled';
% opt.properties = 'constant';    % layers

% Check multiple fluid regions 
opt.grid = 'CR';                  % 'CR' for cappa-rutqvist grid or 'uniform'
opt.inj_rate = 0.002;             % 0.02 for actual CR (unrealistic overpressure)
opt.fluid_model = 'blackoil-ccs';  % 'blackoil-wg', 'blackoil-ccs', 'oil water', 'water'
opt.method = 'fully coupled';
opt.properties = 'layers';   % layers

% Run case
writeIntroText(opt);
optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
optlist = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
[model, states] = runCappaRutqvist2D(optlist{:});


G = createAugmentedGrid(model.G);
figure()
plotNodeDataDeformed(G,sqrt(sum(states{end}.uu.^2,2)),states{end}.uu);
colormap jet;
colorbar;
set(gca, 'YDir','reverse');


figure()
plotToolbar(model.G, states, 'outline', true);
colormap jet;
colorbar;
title(sprintf('fluid model: %s, method: %s', opt.fluid_model, opt.method));
set(gca, 'YDir','reverse');









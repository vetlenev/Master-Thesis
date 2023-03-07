function [G, G_dat, unit, wellId, meshpt] = bMesh(opt, plts, removeS, plotMesh)
%
%
%

%% Get mesh
pth = fullfile(mrstPath('ls-proj'), 'fluidflower/benchmark/mesh');
G_dat  = load(fullfile(pth, opt.mesh{2}));
G = G_dat.G;
wellId = find(~isnan(G_dat.wellNo));
% if removeS
%     unit.ESF = [2, 18, 24];
%     unit.ESFsup = 11;
%     unit.C = [3, 8, 14, 21];
%     unit.D = [5, 9, 20, 22, 30];
%     unit.E = [1, 7, 19, 23, 29, 31];
%     unit.F = [4, 10, 15, 16, 25, 26, 28];
%     unit.G = [6, 12, 13, 17, 27];
% else
%     unit.ESF = [2, 18, 25];
%     unit.ESFsup = 11;
%     unit.C = [3, 8, 14, 22];
%     unit.D = [5, 9, 20, 23, 32];
%     unit.E = [1, 7, 19, 24, 31, 33];
%     unit.F = [4, 10, 15, 16, 27, 28, 30];
%     unit.G = [6, 12, 13, 17, 29];
%     unit.S = [21, 26];
%     
%     % modif G_dat.p for D in fault in original mesh
%     idf = max(G_dat.p) + 1;
%     cf = G_dat.p == 14;
%     cd = all([cf, ...
%               G.cells.centroids(:, 3) < 0.9135, ...
%               G.cells.centroids(:, 3) > 0.714, ...
%               G.cells.centroids(:, 2) > 0.417], 2);
%     G_dat.p(cd) = idf;
%     unit.D = [unit.D idf];
% end
if opt.mesh{1} == 5 || opt.mesh{1} == 10 || opt.mesh{1} == 20 % && v2 only
    if removeS
        unit.ESF = [2, 19, 25];
        unit.ESFsup = 11;
        unit.C = [5, 8, 14, 22];
        unit.D = [4, 9, 17, 21, 23, 31];
        unit.E = [1, 7, 20, 24, 30, 32];
        unit.F = [3, 10, 15, 16, 26, 27, 29];
        unit.G = [6, 12, 13, 18, 28];
    end
else
    if removeS
        unit.ESF = [2, 19, 25];
        unit.ESFsup = 11;
        unit.C = [3, 8, 14, 22];
        unit.D = [5, 9, 17, 21, 23, 31];
        unit.E = [1, 7, 20, 24, 30, 32];
        unit.F = [4, 10, 15, 16, 26, 27, 29];
        unit.G = [6, 12, 13, 18, 28];
    else
        unit.ESF = [2, 19, 26];
        unit.ESFsup = 11;
        unit.C = [3, 8, 14, 23];
        unit.D = [5, 9, 17, 21, 24, 33];
        unit.E = [1, 7, 20, 25, 32, 34];
        unit.F = [4, 10, 15, 16, 28, 29, 31];
        unit.G = [6, 12, 13, 18, 30];
        unit.S = [22, 27];
    end
end


%% meshpt for handy plotting of layer boundaries
meshpt  = load(fullfile(pth, 'benchmark_datapoints_v2.mat'));
meshpt  = meshpt.stratiPoints;


%% Plot?
if plotMesh == 1
   plts.fig3D();
    %%colormap(turbo); plotCellData(G, G_dat.p); view([90, 0]); 
    colr = [0 6 8 15 20 23 26 30];
    colrG = zeros(G.cells.num, 1);
    if ~removeS
        colrG(ismember(G_dat.p, unit.S)) = colr(1);
    end
    colrG(ismember(G_dat.p, unit.ESFsup)) = colr(2);
    colrG(ismember(G_dat.p, unit.ESF)) = colr(3);
    colrG(ismember(G_dat.p, unit.C)) = colr(4);
    colrG(ismember(G_dat.p, unit.D)) = colr(5);
    colrG(ismember(G_dat.p, unit.E)) = colr(6);
    colrG(ismember(G_dat.p, unit.F)) = colr(7);
    colrG(ismember(G_dat.p, unit.G)) = colr(8);
    plotToolbar(G, colrG, 'edgealpha', 1, 'edgecolor', [1 1 1])
    hold on
    plotCellData(G, [1; 1], wellId([1 4]), 'facecolor', 'r');
    plotCellData(G, [1; 1], wellId([2 3]), 'facecolor', [0.8 0.2 0.8]);
    plts.setAxProps(gca), %camlight();
    colormap(copper); %c = colorbar; set(c, 'YTick', sort(colr));
    axis equal off, view([90, 0])
    ax = gca; ax.DataAspectRatio = [0.2 1 1];
    %ylim([0.35 0.55]), zlim([0.35 0.43]) 
end

%% Volume calculation
% Vt = sum(G.cells.volumes)*1e6; %cm3
% V = zeros(numel(uGr), 1);
% for n=1:numel(uGr)
%     % Cell ids of corresponding unit
%     cid = ismember(G_dat.compartID, uGr{n});
%     
%     % volume
%     V(n) = sum(G.cells.volumes(cid));
% end
% V = V*1e6; % cm3

end
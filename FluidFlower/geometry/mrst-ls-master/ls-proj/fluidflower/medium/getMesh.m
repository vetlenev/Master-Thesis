function [G, G_dat, unit, wellId, meshpt] = getMesh(inj_type, opt, plts, plotMesh)
%
%
%

%% Get mesh
if inj_type == 1
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/mesh');
    meshpt  = load(fullfile(pth, 'mesh_datapoints.mat'));
elseif inj_type == 2
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/mesh/AC07');
    meshpt  = load(fullfile(pth, 'mesh_datapoints.mat'));
elseif inj_type == 3
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/mesh/Bilbo');
    meshpt  = load(fullfile(pth, 'bilbo_datapoints.mat'));
end
G_dat  = load(fullfile(pth, opt.mesh{2}));
G = G_dat.G;
wellId = find(~isnan(G_dat.wellNo));
if inj_type == 1
    unit.ESF = [1 16 7];
    unit.C = [3 5 13 15];
    unit.Cf = 11;
    unit.E = [6 8 10 12];
    %unit.Fsup = [2 4 14];
    unit.Fsup = 4;
    unit.Fmid = [2 14];
    unit.Finf = [9 17];       % bottom reservoir
elseif inj_type == 2
    unit.ESF = [1 17 7];
    unit.C = [3 5 14 16];
    unit.Cf = [11 12];
    unit.E = [6 8 10 13];
    %unit.Fsup = [2 4 14];
    unit.Fsup = 4;
    unit.Fmid = [2 15];
    unit.Finf = [9 18];       % bottom reservoir
elseif inj_type == 3
    if opt.mesh{1} == 4 || opt.mesh{1} == 8
        unit.CESF = 8;
        unit.Einf = [6 12];
    elseif opt.mesh{1} == 5 
        unit.CESF = 6;
        unit.Einf = [8 12];
    end
    unit.ESF =[4 9];
    unit.Csup = [1 10]; 
    unit.Cinf = [2 11];
    unit.Esup = 7;
    unit.Finf = 3;
    unit.Fsup = 5;
end

%% meshpt for handy plotting of layer boundaries
meshpt  = meshpt.stratiPoints;

% Depths are positive in G
if inj_type < 3
    meshpt.boundary(:,2) = meshpt.boundary(:,2)*-1;
    for n=1:length(meshpt.lines)
        if n <= numel(meshpt.wells)
            meshpt.wells{n}(:, 2) = meshpt.wells{n}(:, 2)*-1;
        end
        meshpt.lines{n}(:, 2) = meshpt.lines{n}(:, 2)*-1;
    end
end


%% Plot?
if plotMesh == 1
    plts.fig3D();
    Ncompart = max(G_dat.p);
    cmap = copper(Ncompart);
    if inj_type < 3  
        %%colormap(turbo); plotCellData(G, G_dat.p); view([90, 0]); 
        colr = [2, 8, 6, 12, 15, 17];
        colrG = zeros(G.cells.num, 1);
        colrG(ismember(G_dat.p, unit.ESF)) = colr(1);
        colrG(ismember(G_dat.p, unit.C)) = colr(2);
        colrG(ismember(G_dat.p, unit.Cf)) = colr(3);
        colrG(ismember(G_dat.p, unit.E)) = colr(4);
        colrG(ismember(G_dat.p, unit.Fsup)) = colr(5);
        colrG(ismember(G_dat.p, unit.Finf)) = colr(6);
        if isfield(unit, 'Fmid')
            colrG(ismember(G_dat.p, unit.Fmid)) = colr(5);
        end
        plotToolbar(G, colrG, 'edgealpha', 1, 'edgecolor', [0.3 0.3 0.3])
        hold on
        plotCellData(G, [1; 1], wellId([3 5]), 'facecolor', 'r');
        plotCellData(G, [1; 1; 1; 1], wellId([2 4 6 7]), 'facecolor', 'b');
        plotCellData(G, [1; 1], wellId([1 8]), 'facecolor', 'g');
    elseif inj_type == 3
        colr = [2 7 10 12];
        colrG = zeros(G.cells.num, 1);
        colrG(ismember(G_dat.p, unit.ESF)) = colr(1);
        colrG(ismember(G_dat.p, unit.Csup)) = colr(2);
        colrG(ismember(G_dat.p, unit.Cinf)) = colr(2);
        colrG(ismember(G_dat.p, unit.CESF)) = colr(2);
        colrG(ismember(G_dat.p, unit.Esup)) = colr(3);
        colrG(ismember(G_dat.p, unit.Einf)) = colr(3);
        colrG(ismember(G_dat.p, unit.Fsup)) = colr(4);
        colrG(ismember(G_dat.p, unit.Finf)) = colr(4);
        plotToolbar(G, colrG, 'edgealpha', 1, 'edgecolor', [0.3 0.3 0.3])
        hold on
        plotCellData(G, [1; 1], wellId([1 2]), 'facecolor', 'r');
    end
    plts.setAxProps(gca), %camlight();
    colormap(copper); %c = colorbar; set(c, 'YTick', sort(colr));
    axis equal off, view([90, 0])
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
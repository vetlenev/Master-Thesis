%
%
%
% load data
clc, clear, close all force 
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-utils ls-proj
       
%% Read and organize data
% 1. load data for simulation x
dir = 'C:\Users\lsalo\matlab\sim_data\scenario_1\gridRef\mrstdev_gom_adbo_2.5DTriMesh_SC1_CASE1_1mty_spe02perm_45x45x8_pvMult_20yInjTime_20ySimTime\';
%fid = {'ref300'; 'ref200_23LayersRef'; 'ref100_45LayersRef'; 'ref50_91LayersRef'; 'ref25_181LayersRef'};
fid = {'ref300_30LayersRef'; 'ref200_45LayersRef'; 'ref100_90LayersRef'; ...
       'ref50_180LayersRef'; 'ref25_360LayersRef'};
numSim = numel(fid);
r = cell(numSim, 1);
h = waitbar(0, 'Starting to load datasets', 'Name','Loading data...','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
for nSim = 1:numSim
    % Load
    load(strcat(dir, fid{nSim}, '.mat'))
    
    if strcmp(mesh.type, 'ref300')
        resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref300/ucids_resLM2ref.mat');
        elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref300/elemOnTopAreaInX.mat');
    elseif strcmp(mesh.type, 'ref200')
        resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref200/ucids_resLM2ref.mat');
        elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref200/elemOnTopAreaInX.mat');
    elseif strcmp(mesh.type, 'ref100')
        resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref100/ucids_resLM2ref.mat');
        elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref100/elemOnTopAreaInX.mat');
    elseif strcmp(mesh.type, 'ref50')
        resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref50/ucids_resLM2ref.mat');
        elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref50/elemOnTopAreaInX.mat');
    elseif strcmp(mesh.type, 'ref25')
        resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref25/ucids_resLM2ref_upper.mat');
        elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref25/elemOnTopAreaInX_upper.mat');
    end
    cTop = unique(elemOnTop.elemOnTop);
    rLM2 = unique(resLM2ref.resLM2ref);
    
    %% Extract data of interest
    
    % Average gas saturation in refined volume
    ctr                 = max(G.faces.centroids(:,1))/2;
    [minX, maxX]        = deal(ctr-4500, ctr+4500);
    idRef               = all([G.cells.centroids(rLM2, 1) > minX, ...
                               G.cells.centroids(rLM2, 1) < maxX], 2);
    r{nSim}.cRef        = rLM2(idRef);
    r{nSim}.numRef      = numel(rLM2(idRef));
    %r{nSim}.sGavgRef    = mean(states{end}.s(r{nSim}.cRef, 2)); 
    
    % volume(m3) in LM2 with max gas saturation > 0.01
    %idCellSg = states{end}.sMax(r{nSim}.cRef, 2) > 0.01;
    %r{nSim}.volSg    = sum(G.cells.volumes(r{nSim}.cRef(idCellSg)));
    
    % Refined wedge volume
    %r{nSim}.volWedge = sum(G.cells.volumes(r{nSim}.cRef));
    
    % Area with Sg > 0.01
    cTopArea = intersect(cTop, r{nSim}.cRef);                               % cells at the top forming area
    idSg     = states{end}.s(cTopArea, 2) > 0.01;
    f2cn     = gridCellNo(G);                                               % mapping G.cells.faces
    [~, idx] = intersect(f2cn, cTopArea(idSg));                             % find first entry for all cell faces
    numi     = numel(idx);
    j        = repmat([0:4]', numi, 1);
    idx      = repelem(idx, 5) + j;
    fac      = G.cells.faces(idx);                                          % get all faces for top cells in top area
    k        = 5;        facTop = zeros(numi,1);
    for m = 1:numi
        fk = k*m; ik = fk - 4;
        faci = fac(ik:fk);
        [~, idc] = min(G.faces.centroids(faci, 3));
        facTop(m) = faci(idc);                                              % top faces for top cells in top area
    end
    r{nSim}.areaSg = sum(G.faces.areas(facTop));
    
    % Compute value of refined area (should be exactly equal in all cases
    % but it is not)
    clear idx numi j idx fac facTop faci idc
    [~, idx] = intersect(f2cn, cTopArea);
    numi     = numel(idx);
    j        = repmat([0:4]', numi, 1);
    idx      = repelem(idx, 5) + j;
    fac      = G.cells.faces(idx);                                  
    k        = 5;        facTop = zeros(numi,1);
    for m = 1:numi
        fk = k*m; ik = fk - 4;
        faci = fac(ik:fk);
        [~, idc] = min(G.faces.centroids(faci, 3));
        facTop(m) = faci(idc);                                       
    end
    r{nSim}.topArea    = sum(G.faces.areas(facTop));
    r{nSim}.cRefTop    = cTopArea;
    r{nSim}.numRefTop  = numel(cTopArea);
    
    % area (m2) on at the center (y-z plane) of LM2 with MAX gas saturation > 0.01
%     [minX, maxX] = deal(ctr-min(mesh.thick)/2, ctr+min(mesh.thick)/2);
%     idV  = all([G.cells.centroids(rLM2, 1) > minX, ...
%                 G.cells.centroids(rLM2, 1) < maxX], 2);
%     cRefVert = rLM2(idV);                                            % cells at center in ref area
%     idSg     = states{end}.sMax(cRefVert, 2) > 0.01;
%     [~, idx] = intersect(f2cn, cRefVert(idSg));                      % find first entry for all cell faces
%     numi     = numel(idx);
%     j        = repmat([0:4]', numi, 1);
%     idx      = repelem(idx, 5) + j;
%     fac      = G.cells.faces(idx);                                   % get all faces for top cells in top area
%     k        = 5;        facV = zeros(numi,1);
%     for n = 1:numi
%         fk = k*n; ik = fk - 4;
%         faci = fac(ik:fk);
%         [~, idc] = min(G.faces.centroids(faci, 1));
%         facV(n) = faci(idc);                                        % top faces for top cells in top area
%     end
%     r{nSim}.areaSgV = sum(G.faces.areas(unique(facV)));
    
    % Refined area (center)
%     clear idx numi j idx fac facTop faci idc
%     [~, idx] = intersect(f2cn, cRefVert);
%     numi     = numel(idx);
%     j        = repmat([0:4]', numi, 1);
%     idx      = repelem(idx, 5) + j;
%     fac      = G.cells.faces(idx);                                   % get all faces for top cells in top area
%     k        = 5;        facVT = zeros(numi,1);
%     for n = 1:numi
%         fk = k*n; ik = fk - 4;
%         faci = fac(ik:fk);
%         [~, idc] = min(G.faces.centroids(faci, 1));
%         facVT(n) = faci(idc);                                        % top faces for top cells in top area
%     end
%     r{nSim}.sectArea    = sum(G.faces.areas(unique(facVT)));
%     r{nSim}.numRefSect  = numel(facVT);
    
    % plume length
    r{nSim}.plmin = min(G.cells.centroids(r{nSim}.cRef, 2));
    r{nSim}.plmax = max(G.cells.centroids(r{nSim}.cRef, 2));
    
    %% CLEAR NOT NEEDED DATA TO SAVE RAM
    clearvars -except r nSim numSim h fid dir
    
    %% Update waitbar and message
    waitbar(nSim/numSim, h, ['Loaded dataset ' num2str(nSim) ' out of ' num2str(numSim)])
end
delete(h)

%% Plots
cols = jet(nSim);
leg  = fid;
latx = {'Interpreter', 'latex'};
sz   = {'fontSize', 12};
vars = cell2mat(r);
h = [300 200 100 50 25];
%h = [200 100 50 40 25 20 15];

% Areas
subplot(1,2,1)
% plot([vars.numRefTop]/vars(end).topArea, [vars.areaSg], '-sk', ...
%     'linewidth', 1);
plot(1./h, [vars.areaSg], '-sk', 'linewidth', 1);
%legend(leg, latx{:});
%xlabel('$N_\mathrm{c}/$m$^2$', latx{:}, sz{:})
xlabel('$1/h$ [1/m]', latx{:}, sz{:})
ylabel('Area [m$^2$]', latx{:}, sz{:});
grid on
title('Top area with $S_\mathrm{g} > 0.01$', latx{:}, 'fontSize', 11)
ylim([0 1.05*max([vars.areaSg])])

subplot(1,2,2)
plot(1./h, [vars.plmax] - [vars.plmin], '-sk', ...
    'linewidth', 1);
%legend(leg, latx{:});
xlabel('$1/h$ [1/m]', latx{:}, sz{:})
ylabel('Length [m]', latx{:}, sz{:});
grid on
title('Plume length in central section', latx{:}, 'fontSize', 11)
ylim([0 1.05*max([vars.plmax] - [vars.plmin])])

% Disgas
figure(2)
disgas = [2.2879 2.1856 1.9666 1.6484 1.4075]; % Mt; 300, 200, 100, 50, 25
plot(1./h, disgas, '-sk', 'linewidth', 1);
%legend(leg, latx{:});
%xlabel('$N_\mathrm{c}/$m$^2$', latx{:}, sz{:})
xlabel('$1/h$ [1/m]', latx{:}, sz{:})
ylabel('Mass [Mt]', latx{:}, sz{:});
grid on
title('Mass of dissolved gas in the brine', latx{:}, 'fontSize', 11)
ylim([0 1.05*max(disgas)])

% Volumes
% subplot(1,3,3)
% plot([vars.numRef]/vars(end).volWedge, [vars.volSg], '-sk', ...
%     'linewidth', 1);
% %legend(leg, latx{:});
% xlabel('$N_\mathrm{c}/$m$^3$', latx{:}, sz{:})
% ylabel('Volume [m$^3$]', latx{:}, sz{:});
% grid on
% title('Wedge vol. with $S_\mathrm{max,g} > 0.01$', latx{:}, 'fontSize', 11)


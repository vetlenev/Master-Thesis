function plotMeshRes_bo(model, ucids, W, states, state0, ncont, mesh, ...
                        resPlots, opt, injCase)
%
%
%

G = model.G;
rock = model.rock;
fluid = model.fluid;

%% check total CO2 mass
rsSat = fluid.rsSat(states{end}.pressure);
flag = states{end}.rs >= rsSat;
bO = fluid.bO(states{end}.pressure, states{end}.rs, flag);
reg = rock.regions.rocknum;
pvMult = zeros(G.cells.num, 1);
pvMult(reg==1) = fluid.pvMultR{1}(states{end}.pressure(reg==1), state0.pressure(reg==1));
pvMult(reg==2) = fluid.pvMultR{2}(states{end}.pressure(reg==2), state0.pressure(reg==2));
surfVolGas = sum(model.operators.pv .* pvMult .* states{end}.s(:,1) .* states{end}.rs .* bO);
massGasDis = surfVolGas*fluid.rhoGS;
massGasTot = sum(states{end}.FlowProps.ComponentTotalMass{2});
disp(['dissolved mass of CO2 in the brine is ' num2str(massGasDis/10^9) ' Mt'])
disp(['total mass of CO2 in the domain is ' num2str(massGasTot/10^9) ' Mt'])

%% Colors
c_sstsh = [236, 214, 56; 102, 77, 0]./255;
%cmap = flipud(cmap_agu());
cmap = flipud(cmocean('tempo'));
%cmap = hot;
c_sh = [0.5 0.5 0.5];
c_f  = [255, 213, 0]./255;
c_sh_lm2 = [0.5 0.5 0.5];

switch resPlots
    case 'nextToFault'
        %% Reservoir, seal, fault (slice)
        if strcmp(mesh.type, 'coarse') || strcmp(mesh.type, 'ref300') || ...
           strcmp(mesh.type, 'ref200') || strcmp(mesh.type, 'ref100') || ...
           strcmp(mesh.type, 'ref50') || strcmp(mesh.type, 'ref40')
            id = true(G.cells.num, 1);
            id([ucids.unit_cell_ids{1} ucids.unit_cell_ids{2} ucids.unit_cell_ids{3} ...
                ucids.unit_cell_ids{6} ucids.unit_cell_ids{7} ucids.unit_cell_ids{12} ...
                ucids.unit_cell_ids{13}]) = false;
            sealId = 5;
            resId = 4; 
        elseif mesh.reduce == 1
            id = true(G.cells.num, 1);
            id([ucids.unit_cell_ids{3} ucids.unit_cell_ids{4} ucids.unit_cell_ids{9} ...
                ucids.unit_cell_ids{10}]) = false;
            sealId = 2;
            resId = 1;
        elseif mesh.reduce == 2
            id = true(G.cells.num, 1);
            sealId = 2;
            resId = 1;
        end 
        wellLoc = G.cells.centroids(W.cells,:);
        id(G.cells.centroids(:,1) < wellLoc(1)) = false;
        id(G.cells.centroids(:,1) > wellLoc(1)) = false;
        id(G.cells.centroids(:,2) < 1.1*10^4) = false; id(G.cells.centroids(:,2) > 1.7*10^4) = false;
        %id(G.cells.centroids(:,3) < 10^3) = false; id(G.cells.centroids(:,3) > 4*10^3) = false;
        ids = false(G.cells.num, 1);  ids(ucids.unit_cell_ids{2}) = id(ucids.unit_cell_ids{2});
        idf = false(G.cells.num, 1);  idf(ucids.unit_cell_ids{end}) = id(ucids.unit_cell_ids{end});
        figure(90)
        plotToolbar(G, states, ids, 'FaceAlpha', 0, 'EdgeColor', c_sh, ...
            'lineWidth', 0.5, 'EdgeAlpha', 1) % seal edges
        hold on
        plotToolbar(G, states, idf, 'FaceAlpha', 0, 'EdgeColor', c_f, ...
            'lineWidth', 0.5, 'EdgeAlpha', 1) % fault edges
        colormap(cmap)
        plotToolbar(G, states, id, 'FaceAlpha', 0.85);                   % full states plot
        set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
        plotWell(G, W, 'color', 'k', 'height', -1000);
        c = colorbar; caxis([0 1])
        ylim([1.2 1.6]*10^4)
        % outline grid of plotted area can be added with "add patch" button (in GUI)
        
        %% Reservoir top view (plume)
        id2 = false(G.cells.num, 1);
        id2(ucids.unit_cell_ids{resId}) = true;
        %id2(G.cells.centroids(:,2) < 10^4) = false;
        %id2(G.cells.centroids(:,2) > 1.7*10^4) = false;
        figure(91)
        colormap(cmap)
        plotToolbar(G, states, id2, 'FaceAlpha', 1);
        set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
        plotWell(G, W, 'color', 'k', 'height', -1500); c = colorbar; caxis([0 1])
        %ylim([1.2 1.6]*10^4)
        
        %% Fault
%         id3 = false(G.cells.num, 1);
%         id3([ucids.unit_cell_ids{8} ucids.unit_cell_ids{9} ucids.unit_cell_ids{10} ...
%             ucids.unit_cell_ids{11}]) = true;
%         id_lm2_amp = false(G.cells.num, 1); id_lm2_amp(ucids.unit_cell_ids{9}) = true;
%         id_amp = false(G.cells.num, 1); id_amp(ucids.unit_cell_ids{10}) = true;
%         figure(92)
%         plotToolbar(G, states, id_lm2_amp, 'FaceAlpha', 0, 'EdgeColor', c_sh_lm2, ...
%             'lineWidth', 0.5, 'EdgeAlpha', 1) % lm2/amp juxtaposition
%         hold on
%         plotToolbar(G, states, id_amp, 'FaceAlpha', 0, 'EdgeColor', c_sh, ...
%             'lineWidth', 0.5, 'EdgeAlpha', 1) % amp/amp juxtaposition
%         colormap(cmap)
%         plotToolbar(G, states, id3, 'FaceAlpha', 1);
%         set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
%         plotWell(G, W, 'color', 'k', 'height', -1000); c = colorbar; caxis([0 1])
%         ylim([1.2 1.6]*10^4)
        
        %% Fault 2
        cc = max(G.faces.centroids(:,1)/2); limm = 2000;
        idplot1 = G.cells.centroids(ucids.unit_cell_ids{end}, 1)>= cc-limm;
        idplot2 = G.cells.centroids(ucids.unit_cell_ids{end}, 1) <= cc+limm;
        idplot = all([idplot1, idplot2], 2);
        figure(93)
        plotToolbar(G, states, ucids.unit_cell_ids{end}(idplot), 'FaceAlpha', 0.8)
        hold on
        col = ['k','w'];
        for n = 1:numel(ncont)
            idf1 = all([G.nodes.coords(ncont{n}(:,1), 1) >= cc-1.1*limm, ...
                G.nodes.coords(ncont{n}(:,1), 1) <= cc+1.1*limm], 2);
            idf2 = all([G.nodes.coords(ncont{n}(:,2), 1) >= cc-1.1*limm, ...
                G.nodes.coords(ncont{n}(:,2), 1) <= cc+1.1*limm], 2);
            x = G.nodes.coords(ncont{n}(idf1,1), 1); x2 = G.nodes.coords(ncont{n}(idf2,2), 1);
            y = G.nodes.coords(ncont{n}(idf1,1), 2); y2 = G.nodes.coords(ncont{n}(idf2,2), 2);
            z = G.nodes.coords(ncont{n}(idf1,1), 3); z2 = G.nodes.coords(ncont{n}(idf2,2), 3);
            plot3(x, y, z, col(1), 'lineWidth', 5)
            plot3(x2, y2, z2, col(2), 'lineWidth', 5)
        end
        colormap(cmap)
        axis equal off
        view([180, 0])
        
        %% Plot refined area for grid convergence
        if isfield(mesh, 'ref') && mesh.ref == 1
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
            elseif strcmp(mesh.type, 'ref40')
                resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref40/ucids_resLM2ref.mat');
                elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref40/elemOnTopAreaInX.mat');
            elseif strcmp(mesh.type, 'ref25')
                resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref25/ucids_resLM2ref_upper.mat');
                elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref25/elemOnTopAreaInX_upper.mat');
            elseif strcmp(mesh.type, 'ref20')
                resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref20/ucids_resLM2ref_upper.mat');
                elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref20/elemOnTopAreaInX_upper.mat');
            elseif strcmp(mesh.type, 'ref15')
                resLM2ref = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref15/ucids_resLM2ref_upper.mat');
                elemOnTop = load('mrst-dev/mrst-ls/ls-proj/gom/adbo_2.5D/2.5Dmesh/extr_tri/gridRef/ref15/elemOnTopAreaInX_upper.mat');
            end
           rLM2      = unique(resLM2ref.resLM2ref);
           ctr       = max(G.faces.centroids(:,1))/2;
           [minX, maxX] = deal(ctr-4500, ctr+4500);
           idRef  = all([G.cells.centroids(rLM2, 1) > minX, ...
                         G.cells.centroids(rLM2, 1) < maxX], 2);
           cRef      = rLM2(idRef);
           
           figure(94)
           colormap(cmap)
           plotToolbar(G, states, cRef, 'FaceAlpha', 1);
           set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
           plotWell(G, W, 'color', 'k', 'height', -1500); c = colorbar; caxis([0 1])
           
           disp(['Number of elements in refined wedge (HW) is ' num2str(numel(cRef))])
           
           % Refined wedge volume
           volWedge = sum(G.cells.volumes(cRef));
           disp(['Volume of refined wedge is ' ...
                 num2str(volWedge) 'm^3'])
             
           % volume(m3) in LM2 with max gas saturation > 0.01
           idCellSg = states{end}.sMax(cRef, 2) > 0.01;
           volSg    = sum(G.cells.volumes(cRef(idCellSg)));
           disp(['Volume with max Sg > 0.01 in refined wedge is ' ...
                 num2str(volSg) 'm^3'])
           
           % avg gas saturation in refined volume
           sGavgRef = mean(states{end}.s(cRef, 2));
           disp(['Mean gas saturation in refined wedge (HW) after 20y of injection is ' ...
                 num2str(sGavgRef)])
           
           % area (m2) on top of LM2 with gas saturation > 0.01
           cTop     = unique(elemOnTop.elemOnTop);  
           cTopArea = intersect(cTop, cRef);                                % cells at the top forming area
           idSg     = states{end}.s(cTopArea, 2) > 0.01;
           f2cn     = gridCellNo(G);                                        % mapping G.cells.faces
           [~, idx] = intersect(f2cn, cTopArea(idSg));                      % find first entry for all cell faces
           numi     = numel(idx);
           j        = repmat([0:4]', numi, 1);
           idx      = repelem(idx, 5) + j;
           fac      = G.cells.faces(idx);                                   % get all faces for top cells in top area
           k        = 5;        facTop = zeros(numi,1);
           for n = 1:numi
              fk = k*n; ik = fk - 4;
              faci = fac(ik:fk);
              [~, idc] = min(G.faces.centroids(faci, 3));
              facTop(n) = faci(idc);                                        % top faces for top cells in top area
           end
           areaSg = sum(G.faces.areas(facTop));
           disp(['Area with Sg > 0.01 at top of refined wedge is ' ...
                 num2str(areaSg) 'm^2'])
             
           % Refined area
           clear idx numi j idx fac facTop faci idc
           [~, idx] = intersect(f2cn, cTopArea);                         
           numi     = numel(idx);
           j        = repmat([0:4]', numi, 1);
           idx      = repelem(idx, 5) + j;
           fac      = G.cells.faces(idx);                                   % get all faces for top cells in top area
           k        = 5;        facTop = zeros(numi,1);
           for n = 1:numi
              fk = k*n; ik = fk - 4;
              faci = fac(ik:fk);
              [~, idc] = min(G.faces.centroids(faci, 3));
              facTop(n) = faci(idc);                                        % top faces for top cells in top area
           end
           areaSgT = sum(G.faces.areas(facTop));
           disp(['Area of top surface of the refined wedge (full) is ' ...
                 num2str(areaSgT) 'm^2'])       
          
          % area (m2) on at the center (y-z plane) of LM2 with MAX gas saturation > 0.01
           [minX, maxX] = deal(ctr-min(mesh.thick)/2, ctr+min(mesh.thick)/2);
           idV  = all([G.cells.centroids(rLM2, 1) > minX, ...
                         G.cells.centroids(rLM2, 1) < maxX], 2);
           cRefVert = rLM2(idV);                                            % cells at center in ref area
           idSg     = states{end}.sMax(cRefVert, 2) > 0.01;
           [~, idx] = intersect(f2cn, cRefVert(idSg));                      % find first entry for all cell faces
           numi     = numel(idx);
           j        = repmat([0:4]', numi, 1);
           idx      = repelem(idx, 5) + j;
           fac      = G.cells.faces(idx);                                   % get all faces for top cells in top area
           k        = 5;        facV = zeros(numi,1);
           for n = 1:numi
              fk = k*n; ik = fk - 4;
              faci = fac(ik:fk);
              [~, idc] = min(G.faces.centroids(faci, 1));
              facV(n) = faci(idc);                                        % top faces for top cells in top area
           end
           areaSgV = sum(G.faces.areas(facV));
           disp(['Area with max Sg > 0.01 in central section of refined wedge is ' ...
                 num2str(areaSgV) 'm^2'])
           
          % Refined area (center)
           clear idx numi j idx fac facTop faci idc
           [~, idx] = intersect(f2cn, cRefVert);                         
           numi     = numel(idx);
           j        = repmat([0:4]', numi, 1);
           idx      = repelem(idx, 5) + j;
           fac      = G.cells.faces(idx);                                   % get all faces for top cells in top area
           k        = 5;        facTop = zeros(numi,1);
           for n = 1:numi
              fk = k*n; ik = fk - 4;
              faci = fac(ik:fk);
              [~, idc] = min(G.faces.centroids(faci, 1));
              facTop(n) = faci(idc);                                        % top faces for top cells in top area
           end
           areaSgVT = sum(G.faces.areas(facTop));
           disp(['Area of central section of the refined wedge (full) is ' ...
                 num2str(areaSgVT) 'm^2']) 
            
           
        end
        
    case 'nextToFault_sc2'
        %% Reservoir, seal, fault (slice)
        id = false(G.cells.num, 1);
        id([ucids.unit_cell_ids{1:22} ucids.unit_cell_ids{25:32} ucids.unit_cell_ids{58:59}]) = true;
        if opt.sc == 2
            sealId = 2:2:22;
        elseif opt.sc == 1
            sealId = [2:22 58];
        end
        resId = 1;
        
        if strcmp(injCase, 'nextToF')
            lims = [1.2 1.53];
        elseif strcmp(injCase, 'nextToF_FW')
            lims = [1.1 1.4];
        end
        lims2 = [1.25 1.35];
        
        wellLoc = G.cells.centroids(W.cells,:);
        id(G.cells.centroids(:,1) < wellLoc(1)) = false;
        id(G.cells.centroids(:,1) > wellLoc(1)) = false;
        id(G.cells.centroids(:,2) < 9000) = false; id(G.cells.centroids(:,2) > 1.7*10^4) = false;
        %id(G.cells.centroids(:,3) < 10^3) = false; id(G.cells.centroids(:,3) > 4*10^3) = false;
        ids = false(G.cells.num, 1);  ids([ucids.unit_cell_ids{sealId}]) = id([ucids.unit_cell_ids{sealId}]);
        idf = false(G.cells.num, 1);  idf(ucids.unit_cell_ids{end}) = id(ucids.unit_cell_ids{end});
        figure(90)
        plotToolbar(G, states, ids, 'FaceAlpha', 0, 'EdgeColor', c_sh, ...
                    'lineWidth', 0.1, 'EdgeAlpha', 1) % seal edges
        hold on
        plotToolbar(G, states, idf, 'FaceAlpha', 0, 'EdgeColor', 'none') % fault edges
        colormap(cmap)
        plotToolbar(G, states, id, 'FaceAlpha', 0.9);                   % full states plot
        set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
        plotWell(G, W, 'color', 'k', 'height', -1000);
        c = colorbar; caxis([0 1])
        ylim(lims*10^4)
        zlim([1300 2500])
        
        figure(91)
        plotToolbar(G, states, ids, 'FaceAlpha', 0, 'EdgeColor', c_sh, ...
                    'lineWidth', 0.1, 'EdgeAlpha', 1) % seal edges
        hold on
        plotToolbar(G, states, idf, 'FaceAlpha', 0, 'EdgeColor', 'none') % fault edges
        colormap(cmap)
        plotToolbar(G, states, id, 'FaceAlpha', 1);                   % full states plot
        set(gca,'Xdir','reverse'), view(90, 0), camproj perspective; axis equal tight
        plotWell(G, W, 'color', 'k', 'height', -1000);
        c = colorbar; caxis([0 1])
        ylim(lims2*10^4)
        zlim([1600 2200])
        
        % properties on slice
        allVcl = nan(G.cells.num, 1);
        for n = 1:numel(ucids.Vcl)
            if ~isnan(ucids.Vcl(n))
                allVcl(ucids.unit_cell_ids{n}) = ucids.Vcl(n);
            end
        end
        propPlot.kxx_log10 = log10(rock.perm(:, 1)/(milli*darcy));
        propPlot.kyy_log10 = log10(rock.perm(:, 4)/(milli*darcy));
        propPlot.kzz_log10 = log10(rock.perm(:, 6)/(milli*darcy));
        propPlot.poro = rock.poro;
        propPlot.satReg = rock.regions.saturation;
        propPlot.Vcl = allVcl;
        
        f100 = figure(100);
        colormap(copper)
        plotToolbar(G, propPlot, id, 'faceAlpha', 0.8);                   % full states plot
        set(gca,'Xdir','reverse'), axis equal off
        plotWell(G, W, 'color', 'k', 'height', -1000);
        c = colorbar;  c.Label.Interpreter = 'latex'; 
        c.Label.String = '$\log_{10}(k_\mathrm{xx}$ [mD])'; %caxis([0 1])
        c.Label.FontSize = 14;
        ylim(lims*10^4)
        zlim([1300 2500])
        view([90 0]) %view(55, 25)
        set(f100,'position',[50 50 1000 500],...
            'paperunits','points','papersize', [1000 500])
        
        
        %% Reservoir top view (plume)
        id2 = false(G.cells.num, 1);
        id2(ucids.unit_cell_ids{resId}) = true;
        %id2(G.cells.centroids(:,2) < 10^4) = false;
        %id2(G.cells.centroids(:,2) > 1.7*10^4) = false;
        figure(92)
        colormap(cmap)
        plotToolbar(G, states, id2, 'FaceAlpha', 1);
        set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
        plotWell(G, W, 'color', 'k', 'height', -1500); c = colorbar; caxis([0 1])
        
        %% Fault 2
        cc = max(G.faces.centroids(:,1)/2); limm = 2000;
        idplot1 = G.cells.centroids(ucids.unit_cell_ids{61}, 1)>= cc-limm;
        idplot2 = G.cells.centroids(ucids.unit_cell_ids{61}, 1) <= cc+limm;
        idplot = all([idplot1, idplot2], 2);
        figure(93)
        plotToolbar(G, states, ucids.unit_cell_ids{61}(idplot), 'FaceAlpha', 0.8)
        hold on
        col = ['k','w'];
        for n = [1 numel(ncont)-1 numel(ncont)]
            idf1 = all([G.nodes.coords(ncont{n}(:,1), 1) >= cc-1.1*limm, ...
                G.nodes.coords(ncont{n}(:,1), 1) <= cc+1.1*limm], 2);
            idf2 = all([G.nodes.coords(ncont{n}(:,2), 1) >= cc-1.1*limm, ...
                G.nodes.coords(ncont{n}(:,2), 1) <= cc+1.1*limm], 2);
            x = G.nodes.coords(ncont{n}(idf1,1), 1); x2 = G.nodes.coords(ncont{n}(idf2,2), 1);
            y = G.nodes.coords(ncont{n}(idf1,1), 2); y2 = G.nodes.coords(ncont{n}(idf2,2), 2);
            z = G.nodes.coords(ncont{n}(idf1,1), 3); z2 = G.nodes.coords(ncont{n}(idf2,2), 3);
            plot3(x, y, z, col(1), 'lineWidth', 3)
            plot3(x2, y2, z2, col(2), 'lineWidth', 3)
        end
        colormap(cmap)
        axis equal off
        view([150, 15])
        
        
        
    case 'center_meshSc2'
        %% Reservoir, seal, fault (slice)
        id = true(G.cells.num, 1);
        id([ucids.unit_cell_ids{1} ucids.unit_cell_ids{2} ucids.unit_cell_ids{3} ...
            ucids.unit_cell_ids{6} ucids.unit_cell_ids{7} ucids.unit_cell_ids{12} ...
            ucids.unit_cell_ids{13}]) = false;
        wellLoc = G.cells.centroids(W.cells,:);
        id(G.cells.centroids(:,1) < wellLoc(1)) = false;
        id(G.cells.centroids(:,1) > wellLoc(1)) = false;
        id(G.cells.centroids(:,2) < 1.1*10^4) = false; id(G.cells.centroids(:,2) > 2.9*10^4) = false;
        %id(G.cells.centroids(:,3) < 10^3) = false; id(G.cells.centroids(:,3) > 4*10^3) = false;
        ids = false(G.cells.num, 1);  ids(ucids.unit_cell_ids{5}) = id(ucids.unit_cell_ids{5});
        idf = false(G.cells.num, 1);  idf(ucids.unit_cell_ids{end}) = id(ucids.unit_cell_ids{end});
        figure(90)
        plotToolbar(G, states, ids, 'FaceAlpha', 0, 'EdgeColor', c_sh, ...
            'lineWidth', 0.5, 'EdgeAlpha', 1) % seal edges
        hold on
        plotToolbar(G, states, idf, 'FaceAlpha', 0, 'EdgeColor', c_f, ...
            'lineWidth', 0.5, 'EdgeAlpha', 1) % fault edges
        colormap(cmap)
        plotToolbar(G, states, id, 'FaceAlpha', 0.85);                   % full states plot
        set(gca,'Xdir','reverse'), view(55, 12), camproj perspective; axis equal tight
        plotWell(G, W, 'color', 'k', 'height', -1000);
        c = colorbar; caxis([0 1])
        ylim([1.2 2.8]*10^4)
        % outline grid of plotted area can be added with "add patch" button (in GUI)
        
        %% Reservoir top view (plume)
        id2 = false(G.cells.num, 1);
        id2(ucids.unit_cell_ids{4}) = true;
        figure(91)
        colormap(cmap)
        plotToolbar(G, states, id2, 'FaceAlpha', 1);
        set(gca,'Xdir','reverse'), view(55, 25), camproj perspective; axis equal tight
        plotWell(G, W, 'color', 'k', 'height', -1500); c = colorbar; caxis([0 1])
        %ylim([1.2 1.6]*10^4)
        
        %% Fault 2
        cc = max(G.faces.centroids(:,1)/2); limm = 2000;
        idplot1 = G.cells.centroids(ucids.unit_cell_ids{14}, 1)>= cc-limm;
        idplot2 = G.cells.centroids(ucids.unit_cell_ids{14}, 1) <= cc+limm;
        idplot = all([idplot1, idplot2], 2);
        figure(93)
        plotToolbar(G, states, ucids.unit_cell_ids{14}(idplot), 'FaceAlpha', 0.8)
        hold on
        col = ['k','w'];
        for n = 1:numel(ncont)
            idf1 = all([G.nodes.coords(ncont{n}(:,1), 1) >= cc-1.1*limm, ...
                G.nodes.coords(ncont{n}(:,1), 1) <= cc+1.1*limm], 2);
            idf2 = all([G.nodes.coords(ncont{n}(:,2), 1) >= cc-1.1*limm, ...
                G.nodes.coords(ncont{n}(:,2), 1) <= cc+1.1*limm], 2);
            x = G.nodes.coords(ncont{n}(idf1,1), 1); x2 = G.nodes.coords(ncont{n}(idf2,2), 1);
            y = G.nodes.coords(ncont{n}(idf1,1), 2); y2 = G.nodes.coords(ncont{n}(idf2,2), 2);
            z = G.nodes.coords(ncont{n}(idf1,1), 3); z2 = G.nodes.coords(ncont{n}(idf2,2), 3);
            plot3(x, y, z, col(1), 'lineWidth', 5)
            plot3(x2, y2, z2, col(2), 'lineWidth', 5)
        end
        colormap(cmap)
        axis equal off
        view([180, 0])
        
    otherwise
        % do nothing
end

end
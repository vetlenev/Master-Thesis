function plotsSgC(G, rock, model, meshpt, plts, states, sg_bound, c_bound, tid, timesteps)
%
%
%
tsh = cumsum(timesteps)/hour;
plotSgContour = false;
plotSg = true;

% Compute
for n=1:numel(tid)
    sb = states{tid(n)}.s(:,1);
    componentPhaseMass = model.getProp(states{tid(n)}, 'ComponentPhaseMass');
    co2inBrineMass = componentPhaseMass{2,1};  % kg
    Vbrine = G.cells.volumes.*rock.poro.*sb;        % m3
    conc = co2inBrineMass./Vbrine;                  % kg /m3
    sg = states{tid(n)}.s(:,2);
    if n == 1
        Cmax = max(conc);
        Sgmax = max(sg);
    end 

    % CO2 concentration in water
    conc(conc < 1e-9) = 0;
    if plotSgContour
        Fg = scatteredInterpolant(G.cells.centroids(:,2), G.cells.centroids(:,3), sg);
        [L, H] = deal(max(G.nodes.coords(:,2)), max(G.nodes.coords(:,3)));
        [Xq,Yq] = meshgrid(0:.005:L,0:.005:H);
        nx = size(Xq, 2);
        ny = size(Yq, 1);
        Vg = reshape(Fg(Xq(:), Yq(:)), ny, nx);
        h = figure(375);
        Mg = contour(Xq, Yq, Vg, [sg_bound, sg_bound])'; 
        close(h); %set(gca, 'YDir', 'reverse')
        idM = Mg(:,2) <= H;
    end
    
    h = figure(5*n);
    plotCellData(G, conc, 'edgecolor', 'none')
    hold on
    xmx = max(G.faces.centroids(:,1));
    for k=1:numel(meshpt.lines)
        clr = [0.3 0.3 0.3];
        xcrd = repelem(xmx, size(meshpt.lines{k}, 1));
        plot3(xcrd', meshpt.lines{k}(:,1), meshpt.lines{k}(:,2), '-', 'color', clr)
    end
    if plotSgContour
        plot3(xmx*ones(sum(idM), 1), Mg(idM, 1), Mg(idM, 2), '.r', 'markersize', 2)
    end
    %plotLinePath(meshpt.lines(5:end), 'b');
    %plotLinePath(meshpt.lines(1:4), 'r');
    %cmap = getMyCmap('seashore');
    cmap = flipud(cmocean('deep'));
    plts.setAxProps(gca), colormap(cmap), c = colorbar; caxis([0 Cmax])
    axis equal off
    view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
    %set(gca, 'ColorScale', 'log')
    c.Limits = [0 Cmax]; c.Ticks = 0:.1:1.4; 
    ylabel(c, '$C_{\mathrm{CO}_2}$ [kg/m$^3$]', 'fontSize', 14, 'interpreter', 'latex')
    set(h, 'position', [100, 100, 1000, 600])
    if plotSgContour
        exportgraphics(h,['CSg_t' num2str(tsh(tid(n))) '.png'],'ContentType','image',...
                       'Resolution', 300, 'BackgroundColor','w')
    else
        exportgraphics(h,['C_t' num2str(tsh(tid(n))) '.png'],'ContentType','image',...
                       'Resolution', 300, 'BackgroundColor','w')
    end
    close(h);
    
    if plotSg
        h = figure(5*n + 1);
        plotCellData(G, sg, 'edgecolor', 'none')
        hold on
        xmx = max(G.faces.centroids(:,1));
        for k=1:numel(meshpt.lines)
            clr = [0.3 0.3 0.3];
            xcrd = repelem(xmx, size(meshpt.lines{k}, 1));
            plot3(xcrd', meshpt.lines{k}(:,1), meshpt.lines{k}(:,2), '-', 'color', clr)
        end
        cmap = flipud(cmocean('tempo'));
        plts.setAxProps(gca), colormap(cmap), c = colorbar;
        axis equal off
        view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
        set(gca, 'ColorScale', 'log'); caxis([1e-6 1])
        c.Limits = [1e-6 1];
        ylabel(c, '$S_\mathrm{g}$ [-]', 'fontSize', 14, 'interpreter', 'latex')
        set(h, 'position', [100, 100, 1000, 600])
        exportgraphics(h,['Sg_t' num2str(tsh(tid(n))) '.png'],'ContentType','image',...
                       'Resolution', 300, 'BackgroundColor','w')
        close(h);
    end
    
end

end
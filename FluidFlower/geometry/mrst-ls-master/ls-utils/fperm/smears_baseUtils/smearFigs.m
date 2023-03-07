function smearFigs(displayFigures, G, M, f, fw, hw, poroG, permG)
%
%
%

% Plotting utilities
latx = {'Interpreter', 'latex'};
ydim = max(G.faces.centroids(:,2));
rock.poro = poroG;
rock.perm = permG;
idGrid = reshape(transpose(flipud(M.vals)), G.cells.num, 1);

if G.griddim == 2
    if displayFigures(1) == 1
        f1 = figure(1);
        plotToolbar(G, G, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
        axis equal; colorbar; xlim([0 f.T]); ylim([0 f.D]);
        xlabel('x [m]'), ylabel('z [m]');
        title(['Number of cells = ' num2str(G.cells.num)])
        set(f1, 'position', [400, 0, 350, 900]);
    end
    
    if displayFigures(2) == 1
        hh = figure(2);
        h1 = subplot(1,2,1);
        set(h1, 'colormap', jet(M.unit(end)))
        plotToolbar(G, reshape(transpose(flipud(M.units)), G.cells.num, 1), ...
                    'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
        xlim([0 f.T]); ylim([0 f.D]); c = colorbar;
        set(c,'YTick', 1:max(M.unit));
        xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
        title('Unit domains');
        
        
        h2 = subplot(1,2,2);
        layerTops = M.layerTop;
        layerBots = M.layerBot;
        layerCtrs = M.layerCenter;
        cmap = copper(2);
        if unique(idGrid) == 0
            set(h2, 'colormap', cmap(2, :));
        elseif unique(idGrid) == 1
            set(h2, 'colormap', cmap(1, :));
        else
            set(h2, 'colormap', flipud(cmap));
        end
        plotToolbar(G, idGrid, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
        hold on
        for n=1:fw.n
            plot(0,layerTops(n), '+r', 'MarkerSize', 4, 'LineWidth', 1)
            plot(0,layerCtrs(n), 'or', 'MarkerSize', 4, 'LineWidth', 1, 'MarkerFaceColor', 'r')
        end
        plot(0, layerBots(1), '+r', 'MarkerSize', 4, 'LineWidth', 1)
        for n=hw.id(1):hw.id(end)
            plot(f.T, layerTops(n), '+r', 'MarkerSize', 4, 'LineWidth', 1)
            plot(f.T, layerCtrs(n), 'or', 'MarkerSize', 4, 'LineWidth', 1, 'MarkerFaceColor', 'r')
        end
        plot(f.T, layerBots(1), '+r', 'MarkerSize', 4, 'LineWidth', 1)
        hold off
        xlim([0 f.T]); ylim([0 f.D]); c = colorbar; 
        if unique(idGrid) == 0
            caxis([0 0.1]);
        elseif unique(idGrid) == 1
            caxis([0.9 1]);
        end
        %c.Label.Interpreter = 'latex'; c.Label.String = '$\log_{10} k_\mathrm{yy}$ [mD]';
        %c.Label.FontSize = 12;
        xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
        title('Smear location');
        set(hh, 'position', [400, 50, 340, 850]);
        %axis equal tight
    end
    
    if displayFigures(3) == 1
        hh = figure(3);
        h1 = subplot(1,2,1);
        set(h1, 'colormap', copper)
        plotToolbar(G, log10(rock.perm(:,1)/(milli*darcy)), 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
        xlim([0 f.T]); ylim([0 f.D]); c = colorbar;
        caxis([min(log10(rock.perm(:,1)/(milli*darcy))) max(log10(rock.perm(:,1)/(milli*darcy)))]);
        c.Label.Interpreter = 'latex'; %c.Label.String = '$\log_{10} k_\mathrm{xx}$ [mD]';
        c.Label.FontSize = 12;
        xlabel('$x$ [m]', latx{:}); ylabel('$z$ [m]', latx{:})
        title('Across-fault perm'); 
        
        h2 = subplot(1,2,2);
        set(h2, 'colormap', copper)
        plotToolbar(G, log10(rock.perm(:,end)/(milli*darcy)), 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
        xlim([0 f.T]); ylim([0 f.D]); c = colorbar;
        caxis([min(log10(rock.perm(:,1)/(milli*darcy))) max(log10(rock.perm(:,1)/(milli*darcy)))]);
        c.Label.Interpreter = 'latex'; c.Label.String = '$\log_{10} k_\mathrm{zz}$ [mD]';
        c.Label.FontSize = 12;
        set(gca,'fontSize', 13)
        xlabel('$x$ [m]', latx{:}, 'fontSize', 14); ylabel('$z$ [m]', latx{:}, 'fontSize', 14)
        %title('Vertical permeability', 'fontSize', 16);
        set(hh, 'position', [400, 50, 375, 850]);
        %axis equal tight
    end
    
    
elseif G.griddim == 3
    % MRST Grid and grid indexing
    if displayFigures(1) == 1
        f1 = figure(1);
        plotToolbar(G, G, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
        axis equal; colorbar; xlim([0 f.T]); ylim([0 ydim]); zlim([0 f.D]);
        xlabel('x [m]'), ylabel('y [m]'), zlabel('z [m]');
        title(['Number of cells = ' num2str(G.cells.num)])
        set(f1, 'position', [400, 0, 350, 900]); view([-30, 75])
    end
    
    % Smear location
    if displayFigures(2) == 1
        hh = figure(2);
        h1 = subplot(1,2,1);
        set(h1, 'colormap', jet(M.unit(end)))
        %plotToolbar(G, reshape(transpose(flipud(M.units)), G.cells.num, 1), ...
        plotToolbar(G, idGrid.Munits3D, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1);
        xlim([0 f.T]); ylim([0 ydim]); zlim([0 f.D]); c = colorbar;
        set(c,'YTick', 1:max(M.unit));
        xlabel('$x$ [m]', latx{:}); ylabel('$y$ [m]', latx{:}); zlabel('$z$ [m]', latx{:})
        title('Unit domains'); view([-30, 5])
        
        
        h2 = subplot(1,2,2);
        layerTops = f.D-M.layerTop;
        layerBots = f.D-M.layerBot;
        layerCtrs = f.D-M.layerCenter;
        set(h2, 'colormap', flipud(copper(2)))
        plotToolbar(G, idGrid.ids, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
        hold on
        for n=1:fw.n
            plot3([0 0],[0 ydim],[layerTops(n) layerTops(n)], '--r')
            plot3([0 0],[0 ydim],[layerCtrs(n) layerCtrs(n)], '-r', 'lineWidth', 1)
        end
        plot3([0 0],[0 ydim],[layerBots(1) layerBots(1)], '--r')
        for n=hw.id(1):hw.id(end)
            plot3([f.T f.T],[0 ydim],[layerTops(n) layerTops(n)], '--r')
            plot3([f.T f.T],[0 ydim],[layerCtrs(n) layerCtrs(n)], '-r', 'lineWidth', 1)
        end
        plot3([f.T f.T],[0 ydim],[layerBots(1) layerBots(1)], '--r')
        hold off
        xlim([0 f.T]); ylim([0 ydim]); zlim([0 f.D]); c = colorbar;
        %c.Label.Interpreter = 'latex'; c.Label.String = '$\log_{10} k_\mathrm{yy}$ [mD]';
        %c.Label.FontSize = 12;
        xlabel('$x$ [m]', latx{:}); ylabel('$y$ [m]', latx{:}); zlabel('$z$ [m]', latx{:})
        title('Smear location'); view([-30, 5]) %view([-90, 0])
        set(hh, 'position', [400, 50, 600, 850]);
        %axis equal tight
    end
    
    if displayFigures(3) == 1
        hh = figure(3);
        h1 = subplot(1,2,1);
        set(h1, 'colormap', gray(2))
        plotToolbar(G, log10(rock.perm(:,1)/(milli*darcy)), 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
        xlim([0 f.T]); ylim([0 ydim]); zlim([0 f.D]); c = colorbar;
        c.Label.Interpreter = 'latex'; c.Label.String = '$\log_{10} k_\mathrm{xx}$ [mD]';
        c.Label.FontSize = 12;
        xlabel('$x$ [m]', latx{:}); ylabel('$y$ [m]', latx{:}); zlabel('$z$ [m]', latx{:})
        title('Across-fault perm'); view([-30, 5])
        
        h2 = subplot(1,2,2);
        set(h2, 'colormap', gray(2))
        plotToolbar(G, log10(rock.perm(:,6)/(milli*darcy)), 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0);
        xlim([0 f.T]); ylim([0 ydim]); zlim([0 f.D]); c = colorbar;
        c.Label.Interpreter = 'latex'; c.Label.String = '$\log_{10} k_\mathrm{zz}$ [mD]';
        c.Label.FontSize = 12;
        xlabel('$x$ [m]', latx{:}); ylabel('$y$ [m]', latx{:}); zlabel('$z$ [m]', latx{:})
        title('Updip-fault perm'); view([-30, 5])
        set(hh, 'position', [400, 50, 600, 850]);
        %axis equal tight
    end
    
end

end
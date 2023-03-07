function getStatesVideo(dir, type, folderName, idt, plts, ...
                        G, rock, states, model, timesteps, meshpt)
%
%
%

tsm = cumsum(timesteps)'/minute;
v = VideoWriter([dir type '_' folderName '.mp4'], 'MPEG-4');
v.Quality = 100;
v.FrameRate = 1;
open(v);

s = cellfun(@(x) x.s, states(idt), 'UniformOutput', false);
if strcmp(type, 'cgw')
    componentPhaseMass = cellfun(@(x) model.getProp(x, 'ComponentPhaseMass'), ...
                                 states(idt), 'UniformOutput', false);
end

disp('Generating video...')
for k=1:numel(idt)
    
    % Make plot
    plts.fig3D();
    if strcmp(type, 'sg')
        plotCellData(G, s{k}(:,2), 'edgecolor', 'none')
    elseif strcmp(type, 'cgw')
        sb = s{k}(:,1);
        co2inBrineMass = componentPhaseMass{k}{2,1};
        Vbrine = G.cells.volumes.*rock.poro.*sb;
        conc = co2inBrineMass./Vbrine;
        plotCellData(G, conc, 'edgecolor', 'none')
    end
    hold on
    xmx = max(G.faces.centroids(:,1));
    for n=1:numel(meshpt.lines)
        if n < 5
            clr = [0.3 0.3 0.3];
        else
            clr = [0.7 0.7 0.7];
        end
        xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
        plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
    end
    plts.setAxProps(gca)
    %colormap(flipud(cmocean('tempo')))
    colormap(flipud(cmocean('deep')))
    c = colorbar; c.FontSize = 12;
    if strcmp(type, 'sg')
        ylabel(c, '$S_{g}$ [-]', 'fontsize', 18, 'interpreter', 'latex')
    elseif strcmp(type, 'cgw')
        ylabel(c, '$C_{\mathrm{CO}_2}$ [kg/m$^3$]', 'fontSize', 14, 'interpreter', 'latex')
    end
    %caxis([995.1 995.4])
    axis equal off
    view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
    if strcmp(type, 'sg')
         set(gca, 'ColorScale', 'log')
        caxis([1e-6 1])
    elseif strcmp(type, 'cgw')
        caxis([0 1.42])
        %disp(num2str(max(conc)))
    end
    tm = tsm(idt(k));
    h = floor(tm/60);
    m = round(tm - h*60);
    text(max(G.faces.centroids(:,1)), 0.05, 0.05, ...
         [num2str(h) ':' num2str(m) ':00'], ...
         'fontSize', 20, 'fontWeight', 'bold', 'color', 'w');
    
    % Add frame to video
    frame = getframe(gcf);
    writeVideo(v,frame);
    close(gcf)
    
    if mod(k, 10) == 0
        disp([num2str(k) '/' num2str(numel(idt)) ' frames completed.'])
    end
end
close(v)
disp(['Done. Video file saved to ' dir]);
end
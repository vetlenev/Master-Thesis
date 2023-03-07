function plotsBilbo(G, rock, meshpt, plts, s, componentPhaseMass, c_bound, tid, tid2)
%
%
%

% Compute
sb = s{tid}(:,1);
co2inBrineMass = componentPhaseMass{tid}{2,1};  % kg
Vbrine = G.cells.volumes.*rock.poro.*sb;        % m3
conc = co2inBrineMass./Vbrine;                  % kg /m3

sb2 = s{tid2}(:,1);
co2inBrineMass2 = componentPhaseMass{tid2}{2,1};  % kg
Vbrine2 = G.cells.volumes.*rock.poro.*sb2;        % m3
conc2 = co2inBrineMass2./Vbrine2;                  % kg /m3

%% CO2 concentration in water EOI1
conc(conc < 1e-9) = 0;
F = scatteredInterpolant(G.cells.centroids(:,2), G.cells.centroids(:,3), conc);
[L, H] = deal(max(G.nodes.coords(:,2)), max(G.nodes.coords(:,3)));
[Xq,Yq] = meshgrid(0:.005:L,0:.005:H);
nx = size(Xq, 2); 
ny = size(Yq, 1);
Vq = reshape(F(Xq(:), Yq(:)), ny, nx);
h = figure(375);
M = contour(Xq, Yq, Vq, [c_bound, c_bound]); 
close(h); %set(gca, 'YDir', 'reverse')
idM = M(2,:) <= H;

h = figure(1);
plotCellData(G, conc, 'edgecolor', 'none')
hold on
xmx = max(G.faces.centroids(:,1));
for n=1:numel(meshpt.lines)
    clr = [0.3 0.3 0.3];
    xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
    plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
end
%plot3(xmx*ones(sum(idM), 1), M(1,idM), M(2,idM), '-', 'color', 'w')
%plotLinePath(meshpt.lines(5:end), 'b');
%plotLinePath(meshpt.lines(1:4), 'r');
%cmap = getMyCmap('seashore');
cmap = flipud(cmocean('deep'));
plts.setAxProps(gca), colormap(jet), c = colorbar; %clim([0 40000])
axis equal off
view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
%set(gca, 'ColorScale', 'log')
c.Limits = [0 1.4]; c.Ticks = 0:.1:1.4; 
ylabel(c, '$C_{\mathrm{CO}_2}$ [kg/m$^3$]', 'fontSize', 14, 'interpreter', 'latex')
set(h, 'position', [100, 100, 1000, 600])

%% CO2 concentration in water EOI2
conc2(conc2 < 1e-9) = 0;
F2 = scatteredInterpolant(G.cells.centroids(:,2), G.cells.centroids(:,3), conc2);
Vq2 = reshape(F2(Xq(:), Yq(:)), ny, nx);
h = figure(375);
M2 = contour(Xq, Yq, Vq2, [c_bound, c_bound]); 
close(h); %set(gca, 'YDir', 'reverse')
idM2top = M2(2,:) <= 0.325;
idM2bot = all([M2(2,:)' > 0.325, M2(2,:)' <= H], 2); 

h = figure(2);
plotCellData(G, conc2, 'edgecolor', 'none')
hold on
xmx = max(G.faces.centroids(:,1));
for n=1:numel(meshpt.lines)
    clr = [0.3 0.3 0.3];
    xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
    plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
end
%plot3(xmx*ones(sum(idM2top), 1), M2(1,idM2top), M2(2,idM2top), '-', 'color', 'w')
%plot3(xmx*ones(sum(idM2bot), 1), M2(1,idM2bot), M2(2,idM2bot), '-', 'color', 'w')
%plotLinePath(meshpt.lines(5:end), 'b');
%plotLinePath(meshpt.lines(1:4), 'r');
%cmap = getMyCmap('seashore');
cmap = flipud(cmocean('deep'));
plts.setAxProps(gca), colormap(jet), c = colorbar; %clim([0 40000])
axis equal off
view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
%set(gca, 'ColorScale', 'log')
c.Limits = [0 1.4]; c.Ticks = 0:.1:1.4; 
ylabel(c, '$C_{\mathrm{CO}_2}$ [kg/m$^3$]', 'fontSize', 14, 'interpreter', 'latex')
set(h, 'position', [100, 100, 1000, 600])

%% Gas saturation EOI1
h = figure(3);
plotCellData(G, s{tid}(:,2), 'edgecolor', 'none')
hold on
xmx = max(G.faces.centroids(:,1));
for n=1:numel(meshpt.lines)
    clr = [0.3 0.3 0.3];
    xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
    plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
end
cmap = flipud(cmocean('tempo'));
plts.setAxProps(gca), colormap(cmap), c = colorbar; %clim([0 40000])
axis equal off
view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
set(gca, 'ColorScale', 'log')
c.Limits = [1e-6 1];
ylabel(c, '$S_\mathrm{g}$ [-]', 'fontSize', 14, 'interpreter', 'latex')
set(h, 'position', [100, 100, 1000, 600])

h = figure(4);
plotCellData(G, s{tid2}(:,2), 'edgecolor', 'none')
hold on
xmx = max(G.faces.centroids(:,1));
for n=1:numel(meshpt.lines)
    clr = [0.3 0.3 0.3];
    xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
    plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
end
cmap = flipud(cmocean('tempo'));
plts.setAxProps(gca), colormap(cmap), c = colorbar; %clim([0 40000])
axis equal off
view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
set(gca, 'ColorScale', 'log')
c.Limits = [1e-6 1];
ylabel(c, '$S_\mathrm{g}$ [-]', 'fontSize', 14, 'interpreter', 'latex')
set(h, 'position', [100, 100, 1000, 600])


end
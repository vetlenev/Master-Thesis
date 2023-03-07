function [bc, model] = bBC(G, fluid, model, p_r, opt)
%
%
%
g = norm(gravity);
if strcmp(opt.topboundary, 'straight')
    f = false(G.faces.num, 1);
    f(boundaryFaces(G)) = true;
    f = all([G.faces.centroids(:, 3) < 1e-4, f], 2);
    if strcmp(opt.bctype, 'cp')
        bc = addBC([], find(f), 'pressure', p_r, 'sat', [1 0]);
        %     fig3D(); plotFaces(G); plotFaces(G, f, 'edgecolor', 'r')
        %     setAxProps(gca), colormap(jet), c = colorbar;
        %     axis equal off
    elseif strcmp(opt.bctype, 'pvm')    
        cellsext = unique(reshape(G.faces.neighbors(f, :), [], 1));
        cellsext(cellsext==0) = [];
        %     fig3D(); plotGrid(G); hold on, plotGrid(G, cellsext, 'edgecolor', 'r')
        %     setAxProps(gca), axis equal off, view([90 0]);
        model.operators.pv(cellsext) = model.operators.pv(cellsext)*1e4;
        bc = [];
    end
else
    % Find top faces
    f = false(G.faces.num, 1);
    f(boundaryFaces(G)) = true;
    zmean = mean(G.cells.centroids(:, 3));
    f = all([G.faces.centroids(:, 3) < zmean, f], 2);
    f = all([f, abs(G.faces.normals(:, 3)) > abs(0.9*G.faces.normals(:, 2))], 2);
    f = find(all([f, abs(G.faces.normals(:, 3)) > abs(0.9*G.faces.normals(:, 1))], 2));
    f(G.faces.centroids(f,3)>0.1) = [];
    
    
%     plts.fig3D(); plotFaces(G); plotFaces(G, f, 'edgecolor', 'r', 'facecolor', 'none')
%     plts.setAxProps(gca), colormap(jet), c = colorbar;
%     axis equal off

    if strcmp(opt.bctype, 'cp')
        [z_0, z_max] = deal(min(G.faces.centroids(f,3)), max(G.faces.centroids(f,3)));
        equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [z_0, z_max], p_r);
        fp_val = reshape(deval(equil, G.faces.centroids(f,3)), [], 1);  clear equil
        bc = addBC([], f, 'pressure', fp_val, 'sat', [1, 0]);
        
        %cells_top = unique([G.faces.neighbors(f,1); G.faces.neighbors(f,2)]);
        %cells_top(cells_top==0) = [];
        %model.operators.pv(cells_top) = model.operators.pv(cells_top)*1e3;
        
        %plot3(G.faces.centroids(f,1), G.faces.centroids(f,2), (1-G.faces.centroids(f,3)/3), '-k')
        %hold on
        %plot3(G.faces.centroids(f,1), G.faces.centroids(f,2), fp_val/max(fp_val), '-r')
    elseif strcmp(opt.bctype, 'pvm')
        cellsext = unique(reshape(G.faces.neighbors(f, :), [], 1));
        cellsext(cellsext==0) = [];
        %plts.fig3D(); plotGrid(G); hold on, plotGrid(G, cellsext, 'edgecolor', 'r')
        %plts.setAxProps(gca), axis equal off, view([90 0])
        model.operators.pv(cellsext) = model.operators.pv(cellsext)*10^9; 
        bc = [];
    end
end

end
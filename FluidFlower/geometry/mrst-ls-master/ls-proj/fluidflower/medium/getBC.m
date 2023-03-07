function [bc] = getBC(G, fluid, p_r, opt)
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
    f = all([G.faces.centroids(:, 3) < 0.05, f], 2);
    f = all([f, abs(G.faces.normals(:, 3)) > abs(0.9*G.faces.normals(:, 2))], 2);
    f = find(all([f, abs(G.faces.normals(:, 3)) > abs(0.9*G.faces.normals(:, 1))], 2));

%     plts.fig3D(); plotFaces(G); plotFaces(G, f, 'edgecolor', 'r')
%     plts.setAxProps(gca), colormap(jet), c = colorbar;
%     axis equal off; view([90 0])

    if strcmp(opt.bctype, 'cp')
        [z_0, z_max] = deal(min(G.faces.centroids(f,3)), max(G.faces.centroids(f,3)));
        equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [z_0, z_max], p_r);
        fp_val = reshape(deval(equil, G.faces.centroids(f,3)), [], 1);  clear equil
        bc = addBC([], f, 'pressure', fp_val, 'sat', [1, 0]);
    elseif strcmp(opt.bctype, 'pvm')
        cellsext = unique(reshape(G.faces.neighbors(f, :), [], 1));
        cellsext(cellsext==0) = [];
        %fig3D(); plotGrid(G); hold on, plotGrid(G, cellsext, 'edgecolor', 'r')
        %setAxProps(gca), axis equal off, view([90 0])
        model.operators.pv(cellsext) = model.operators.pv(cellsext)*1e4; 
        bc = [];
    end
end

end
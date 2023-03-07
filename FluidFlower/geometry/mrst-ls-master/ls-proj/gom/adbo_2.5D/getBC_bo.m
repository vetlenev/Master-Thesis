function [bc, model] = getBC_bo(G, fluid, model, opt, mesh, ucids)
%
%
%

switch opt.bc
    case 'mixed'
        g = norm(gravity);
        rho_wr = fluid.rhoOS*kilogram/meter^3;
        p_r = g*rho_wr*opt.watCol;
        f = boundaryFaces(G);
        fp = f(G.faces.centroids(f,3) < max(G.faces.centroids(f,3))-5);
        fp =fp(G.faces.centroids(fp,3) > min(G.faces.centroids(fp,3))+5);
        [z_0, z_max] = deal(min(G.faces.centroids(fp,3)), max(G.faces.centroids(fp,3)));
        equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS(1), [z_0, z_max], p_r);
        fp_val = reshape(deval(equil, G.faces.centroids(fp,3)), [], 1);  clear equil
        bc = addBC([], fp, 'pressure', fp_val, 'sat', [1, 0]);
        bc.dissolution = [1, 0; 0, 0];  % H2O+Salt and CO2 fractions (cols.) in aqueous and gas phases (row of matrix)
    case 'noflow'
        bc = [];
        
    case 'pvMult' 
        % this is applied to the cells in the perimeter of the
        % injection layer
        if mesh.reduce == 0
            id_resLM2 = 4;
        else
            id_resLM2 = 1;
        end
        bc = [];
        %idext1 = G.faces.neighbors(:,2) == 0; idext2 = G.faces.neighbors(:,1) == 0;
        %cellsext = [G.faces.neighbors(idext1, 1); G.faces.neighbors(idext2, 2)];
        lim = [max(mesh.thick) max(G.faces.centroids(:,2))-max(mesh.thick);
               150 44900]; 
        % note that lim for y dim would not work for a mesh with different
        % triangle size (other than 100m) towards the reservoir extremes.
        cells = sort(ucids.unit_cell_ids{id_resLM2});
        idx = [G.cells.centroids(cells, 1) < lim(1,1) ...
               G.cells.centroids(cells, 1) > lim(1,2)];
        idy = [G.cells.centroids(cells, 2) < lim(2,1) ...
               G.cells.centroids(cells, 2) > lim(2,2)];
        ids = any([any(idx, 2) any(idy, 2)], 2);
        cellsext = cells(ids);
        model.operators.pv(cellsext) = model.operators.pv(cellsext)*1e4;
        
    otherwise
        error('This type of boundary condition is not supported')
end

% % Plot
% clf, plotGrid(G, cells, 'FaceColor', 'none', 'EdgeColor', [0.8 0.8 0.8]); set(gca,'Xdir','reverse'); 
% view([55,25]), camproj perspective; axis equal tight;
% plotGrid(G, cellsext, 'FaceColor', 'r', 'EdgeColor', [0.8 0.8 0.8]);
% 
% plotFaces(G, fp, fp_val/barsa, 'EdgeColor', 'r'); colorbar;

end
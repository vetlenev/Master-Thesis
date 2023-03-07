function [areas_A, areas_B] = getAreas(G, model, timesteps, states, problem, saveCsv)

ts = cumsum(timesteps);
rock = model.rock;
% Get interpolant from benchmarkMesh, for thickness (should have saved it)
% You can also get thickness from nodes in G structure but more tedious to
% then compute the cell thickness.
x = [0 0.03 0.73 1.33 2.13 2.83 2.86];
ymax = max(G.nodes.coords(:,3));
y = [ymax-0.03 ymax-0.33 ymax-0.63 ymax-0.93 ymax-1.23 ymax-1.53];
[X, Y] = meshgrid(x, y);
v = [0.019 0.019 0.02 0.019 0.02 0.019 0.019;
    0.019 0.019 0.024 0.024 0.024 0.019 0.019;
    0.019 0.019 0.026 0.028 0.026 0.019 0.019;
    0.019 0.019 0.027 0.025 0.026 0.019 0.019;
    0.019 0.019 0.026 0.026 0.024 0.019 0.019;
    0.019 0.019 0.019 0.023 0.02 0.019 0.019];
h = 5/1000; % 5mm mesh
[Xq,Yq] = meshgrid(0:h:2.86,0:h:ymax);
Vq = interp2(X,Y,v,Xq,Yq,'spline');
%surf(Xq,Yq,Vq,'edgecolor','none'); view([0 -90]);
F = scatteredInterpolant(Xq(:),Yq(:),Vq(:));
    
% Get other data                    
lim_A = [1.13 2.83; 0.71 1.31];   % [y; z]
id_A = all([G.cells.centroids(:,2) >= lim_A(1) ...
            G.cells.centroids(:,2) <= lim_A(1,2) ...
            G.cells.centroids(:,3) >= lim_A(2,1), ...
            G.cells.centroids(:,3) <= lim_A(2,2)], 2);
area_A_cm2 = diff(lim_A(1,:))*diff(lim_A(2,:))*1e4;    % cm2
lim_B = [0.03 1.13; 0.11 0.71];   % [y; z]
id_B = all([G.cells.centroids(:,2) >= lim_B(1) ...
            G.cells.centroids(:,2) <= lim_B(1,2) ...
            G.cells.centroids(:,3) >= lim_B(2,1), ...
            G.cells.centroids(:,3) <= lim_B(2,2)], 2);
area_B_cm2 = diff(lim_B(1,:))*diff(lim_B(2,:))*1e4;
thick = F(G.cells.centroids(:,2), G.cells.centroids(:,3));
sg_bound = 1e-3;
c_bound = 0.15*1.4;
areas_A = nan(numel(states), 3);
areas_B = nan(numel(states), 3);
for n=1:numel(states)
    sb = states{n}.s(:,1);
    componentPhaseMass = model.getProp(states{n}, 'ComponentPhaseMass');
    co2inBrineMass = componentPhaseMass{2,1};       % kg
    Vbrine = G.cells.volumes.*rock.poro.*sb;        % m3
    conc = co2inBrineMass./Vbrine;                  % kg /m3
    idAboveG = all([states{n}.s(:,2) > sg_bound, id_A], 2);
    A_with_co2gA_cm2 = sum(G.cells.volumes(idAboveG) ./ thick(idAboveG))*1e4;  % Sum of cell areas
    idAboveC = all([conc > c_bound, id_A], 2);
    A_with_co2bA_cm2 = sum(G.cells.volumes(idAboveC) ./ thick(idAboveC))*1e4;
    A_with_watA_cm2 = area_A_cm2 - A_with_co2bA_cm2;    % A_with_co2b includes Agas
    areas_A(n,:) = [A_with_watA_cm2 A_with_co2bA_cm2 A_with_co2gA_cm2];
    
    idAboveG = all([states{n}.s(:,2) > sg_bound, id_B], 2);
    A_with_co2gB_cm2 = sum(G.cells.volumes(idAboveG) ./ thick(idAboveG))*1e4;  % Sum of cell areas
    idAboveC = all([conc > c_bound, id_B], 2);
    A_with_co2bB_cm2 = sum(G.cells.volumes(idAboveC) ./thick(idAboveC))*1e4;
    A_with_watB_cm2 = area_B_cm2 - A_with_co2bB_cm2;
    areas_B(n,:) = [A_with_watB_cm2 A_with_co2bB_cm2 A_with_co2gB_cm2];
end

% write table
if saveCsv
    hdrs = {'t_s', 'A_w_cm2', 'A_co2b_cm2', 'A_co2g_cm2', 'B_w_cm2', 'B_co2b_cm2', 'B_co2g_cm2',};
    vartyp  = cellstr(repmat('double', numel(hdrs), 1));
    nr = numel(states);
    t = table('Size', [nr, numel(hdrs)], 'VariableTypes', vartyp, 'VariableNames', hdrs);
    t.t_s = ts';
    t.A_w_cm2 = areas_A(:,1);
    t.A_co2b_cm2 = areas_A(:,2);
    t.A_co2g_cm2 = areas_A(:,3);
    t.B_w_cm2 = areas_B(:,1);
    t.B_co2b_cm2 = areas_B(:,2);
    t.B_co2g_cm2 = areas_B(:,3);
    writetable(t, ['areas_' problem.Name '.csv'], 'Delimiter', ',');
end
    
end
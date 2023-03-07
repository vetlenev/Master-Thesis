function t_fingers_s = getFingerTimes(G, G_dat, model, states, unit, timesteps, ...
                                      c_bound, name, saveCsv)
%
%
%
rock = model.rock;
ts = cumsum(timesteps);
yrang = [0.05 0.2; % bottom
         0.1 0.4];  % mid C
H = max(G.cells.centroids(:,3));

id_cells_bot = all([G.cells.centroids(:,3) > (H-0.002), ...
                    G.cells.centroids(:,2) >= yrang(1,1), ...
                    G.cells.centroids(:,2) <= yrang(1,2)], 2);
%plotGrid(G, 'facecolor', 'none'); hold on; plotGrid(G, id_cells_bot); view([90 0])

id_cells_mid = all([G.cells.centroids(:,2) >= yrang(2,1), ...
                    G.cells.centroids(:,2) <= yrang(2,2), ...
                    G_dat.p==unit.C(2)], 2);
%plotGrid(G, 'facecolor', 'none'); hold on; plotGrid(G, id_cells_mid); view([90 0])

t_fingers_s = nan(2,1);
check_bot = true;
check_mid = true;
for n=150:numel(states)
    sb = states{n}.s(:,1);
    componentPhaseMass = model.getProp(states{n}, 'ComponentPhaseMass');
    co2inBrineMass = componentPhaseMass{2,1};   % kg
    Vbrine = G.cells.volumes.*rock.poro.*sb;    % m3
    conc = co2inBrineMass./Vbrine;              % kg /m3
    if any(conc(id_cells_bot) > c_bound(1)) && check_bot
        t_fingers_s(1) = ts(n);
        check_bot = false;
    end
    if any(conc(id_cells_mid) > c_bound(2)) && check_mid
        t_fingers_s(2) = ts(n);
        check_mid = false;
    end
    if ~check_bot && ~check_mid
        break
    end
end
%assert(t_fingers_s(2) > t_fingers_s(1))

if saveCsv
    hdrs = {'measure', 't_s'};
    vartyp  = cellstr(['string'; repmat('double', numel(hdrs)-1, 1)]);
    nr = 2;
    t = table('Size', [nr, numel(hdrs)], 'VariableTypes', vartyp, 'VariableNames', hdrs);
    t.measure = {'t_finger_bot'; 't_finger_mid'};
    t.t_s = {t_fingers_s(1); t_fingers_s(2)};
    writetable(t, ['tFingers_' name '.csv'], 'Delimiter', ',');
end

end
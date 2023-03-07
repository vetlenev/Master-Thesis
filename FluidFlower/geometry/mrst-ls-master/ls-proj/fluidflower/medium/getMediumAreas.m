function A_simulation_cm2 = getMediumAreas(G, G_dat, unit, inj_type, states, model, ...
                                           timesteps, sg_bound, c_bound, name, saveCsv)
%
%
%
ts = cumsum(timesteps);
rock = model.rock;
if inj_type == 1
    tid = 54;   %24, 54
elseif inj_type == 2
    tid = 290;
elseif inj_type == 3
   tid = 70;
end
idcF = G_dat.p == unit.Finf(1);
sb = states{tid}.s(:,1);
componentPhaseMass = model.getProp(states{tid}, 'ComponentPhaseMass');
co2inBrineMass = componentPhaseMass{2,1};  % kg
Vbrine = G.cells.volumes.*rock.poro.*sb;        % m3
conc = co2inBrineMass./Vbrine;                  % kg /m3
V_with_co2g = sum(G.cells.volumes(states{tid}.s(idcF,2) >= sg_bound)); % Sum of cell volumes
%V_with_co2b = sum(G.cells.volumes(states{tid}.PVTProps.Density{1} > 995.18));
V_with_co2b = sum(G.cells.volumes(all([conc > c_bound, idcF], 2)));
%Vco2_g_cells = model.operators.pv .* states{tid}.s(:,2);
%Vco2_g = sum(Vco2_g_cells);                                     % Actual g volume
%Mco2_g = sum(Vco2_g_cells.*states{tid}.PVTProps.Density{2});
%Mco2_tot = sum(states{tid}.FlowProps.ComponentTotalMass{2});
%Mco2_b = Mco2_tot - Mco2_g;
%rsSat = fluid.rsSat(states{tid}.pressure);
%flag = states{tid}.rs >= rsSat;
%bO = fluid.bO(states{tid}.pressure, states{tid}.rs, flag);
%Vco2_b = sum(model.operators.pv .* states{tid}.s(:,1) .* states{tid}.rs .* bO);
%Vco2_tot = Vco2_b + Vco2_g;
%Mco2_b_check = Vco2_b*fluid.rhoGS;
%disp(['Surface V of dissolved CO2 is ' num2str(Vco2_b*10^6) ' ml'])
%disp(['Mass of CO2 in the brine is ' num2str(Mco2_b*10^6) ' mg'])
%disp(['total V of CO2 in the domain is ' num2str(Vco2_tot*10^6) ' ml'])
%disp(['total mass of CO2 in the domain is ' num2str(Mco2_tot*10^6) ' mg'])

% End of I2
if inj_type == 1
    tid2 = 154; %67, 154;
    idc = G_dat.p==4;
    idc2 = G_dat.p == 2;
elseif inj_type == 2
    tid2 = 290;
    idc = G_dat.p==4;
    idc2 = G_dat.p == 2;
    idc3 = G_dat.p == 15;
elseif inj_type == 3 
    tid2 = 153;
    idc = G_dat.p == unit.Fsup;
    idc2 = G_dat.p == unit.Esup;
end
sb = states{tid2}.s(:,1);
componentPhaseMass = model.getProp(states{tid2}, 'ComponentPhaseMass');
co2inBrineMass = componentPhaseMass{2,1}; % kg
Vbrine = G.cells.volumes.*rock.poro.*sb;        % m3
conc = co2inBrineMass./Vbrine;                  % kg /m3
V_with_co2g2 = sum(G.cells.volumes(states{tid2}.s(idc,2) >= sg_bound)); % Sum of cell volumes
% V_with_co2b2 = sum(G.cells.volumes(all([states{tid}.PVTProps.Density{1} > ...
%                                        995.15, G_dat.p==4], 2)));
V_with_co2b2 = sum(G.cells.volumes(all([conc > c_bound, idc], 2)));

V_with_co2g_mid = sum(G.cells.volumes(states{tid2}.s(idc2,2) >= 1e-3)); % Sum of cell volumes
% V_with_co2b_mid = sum(G.cells.volumes(all([states{tid2}.PVTProps.Density{1} > ...
%                                        995.15, G_dat.p==2], 2)));
V_with_co2b_mid = sum(G.cells.volumes(all([conc > c_bound, idc2], 2)));
%plotGrid(G, 'facecolor', 'none'); plotGrid(G, all([conc > c_bound, idcF], 2)); view([90 0])
if inj_type == 2
    V_with_co2g_mid2 = sum(G.cells.volumes(states{tid2}.s(idc3,2) >= 1e-3));
    V_with_co2b_mid2 = sum(G.cells.volumes(all([conc > c_bound, idc3], 2)));
end

% Areas
if inj_type == 1 || inj_type == 3
    V_simulation = [V_with_co2g V_with_co2g_mid V_with_co2g2 ...
                    V_with_co2b V_with_co2b_mid V_with_co2b2];
elseif inj_type == 2
    V_simulation = [V_with_co2g V_with_co2g_mid V_with_co2g_mid2 V_with_co2g2 ...
                    V_with_co2b V_with_co2b_mid V_with_co2b_mid2 V_with_co2b2];
end
w = max(G.nodes.coords(:,1));
A_simulation_cm2 = V_simulation/w*1e4;

% write table
if saveCsv
    n1 = ['Area_cm2_t' num2str(ts(tid))];
    n2 = ['Area_cm2_t' num2str(ts(tid2))];
    if inj_type == 2
        n2 = ['Area_cm2_t' num2str(ts(tid2)) '_'];
    end
    hdrs = {'A_name', n1, n2};
    vartyp  = cellstr(['string'; repmat('double', numel(hdrs)-1, 1)]);
    nr = numel(A_simulation_cm2);
    t = table('Size', [nr, numel(hdrs)], 'VariableTypes', vartyp, 'VariableNames', hdrs);
    if inj_type == 1 || inj_type == 3
        t.A_name = {'V_with_co2g'; 'V_with_co2g_mid'; 'V_with_co2g2'; ...
                    'V_with_co2b'; 'V_with_co2b_mid'; 'V_with_co2b2'};
        t(:,2) = num2cell([A_simulation_cm2(1); nan(2, 1); ...
                           A_simulation_cm2(4); nan(2,1)]);
        t(:,3) = num2cell([nan; A_simulation_cm2(2:3)'; ...
                           nan; A_simulation_cm2(5:6)']);
    elseif inj_type == 2
        t.A_name = {'V_with_co2g'; 'V_with_co2g_mid'; 'V_with_co2g_mid2'; 'V_with_co2g2'; ...
                    'V_with_co2b'; 'V_with_co2b_mid'; 'V_with_co2b_mid2'; 'V_with_co2b2'};
        t(:,2) = num2cell([A_simulation_cm2(1); nan(3, 1); ...
                           A_simulation_cm2(5); nan(3,1)]);
        t(:,3) = num2cell([nan; A_simulation_cm2(2:4)'; nan; A_simulation_cm2(6:8)']);
    end
    writetable(t, ['areas_' name '.csv'], 'Delimiter', ',');
end

end
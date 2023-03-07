function [total_mass] = getMass(model, states, timesteps, name, saveCsv)
%
%
%

ts = cumsum(timesteps);
total_mass = nan(numel(states), 4);
for n=1:numel(states)
    componentPhaseMass = model.getProp(states{n}, 'ComponentPhaseMass');
    relativePermeability = model.getProp(states{n}, 'RelativePermeability');
    componentTotalMass = model.getProp(states{n}, 'ComponentTotalMass');
    id_mob = relativePermeability{2} >= 1e-9;
    mob = sum(componentPhaseMass{2,2}(id_mob));
    immob = sum(componentPhaseMass{2,2}(~id_mob));
    diss = sum(componentPhaseMass{2,1});
    tot = sum(componentTotalMass{2}); 
    assert(abs(tot-(mob+immob+diss)) < 1e-9)
    
    total_mass(n,:) = [mob, immob, diss, tot];
end

if saveCsv
    hdrs = {'t_s', 'mob_kg', 'immob_kg', 'diss_kg', 'tot_kg'};
    vartyp  = cellstr(repmat('double', numel(hdrs), 1));
    nr = numel(states);
    t = table('Size', [nr, numel(hdrs)], 'VariableTypes', vartyp, ...
              'VariableNames', hdrs);
    t(:,1) = num2cell(ts');
    t(:,2) = num2cell(total_mass(:,1));
    t(:,3) = num2cell(total_mass(:,2));
    t(:,4) = num2cell(total_mass(:,3));
    t(:,5) = num2cell(total_mass(:,4));
    writetable(t, ['mass_' name '.csv'], 'Delimiter', ',');
end




end
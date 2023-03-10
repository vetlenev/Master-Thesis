% Assume that this repo is cloned in the same spot as the official
% fluidsolver/spe11 repo and that we are currently in spe11-utils/utils
if false
    spe11_folder = fullfile(pwd(), '..', '..', 'spe11');
else
    spe11_folder = fullfile(pwd(), 'Master-Thesis', 'FluidFlower', 'geometry', 'deck');
end
assert(isfolder(spe11_folder));
%%
do_plot = false;
np = 50;
p_min = 1*atm;
p_max = 2.5*barsa;
T = 20 + 273.15;
nusat = 5;

dissolve_gas = true;
vaporize_water = false;

clc
fprintf('To generate, run the following:\n');
for i = 1:2
    if i == 1
        name = '-cH2O';
    else
        name = '-cCO2';
    end
    fprintf('python ./make_component_table.py -t1 %f -t2 %f -nt 2 -p1 %f -p2 %f -np %d %s\n', T, 1.1*T, p_min, p_max, np, name)
end
fprintf('python ./make_solubility_table.py -t1 %f -t2 %f -nt 2 -p1 %f -p2 %f -np %d\n', T, 1.1*T, p_min, p_max, np)

tab_h2o = readtable(fullfile(spe11_folder, 'thermodynamics', 'h2ovalues.csv'));
tab_co2 = readtable(fullfile(spe11_folder, 'thermodynamics', 'co2values.csv'));
tab_sol = readtable(fullfile(spe11_folder, 'thermodynamics', 'solubilities.csv'));


writeFluidFlowerPROPS(tab_h2o, tab_co2, tab_sol, 'rs', true, 'rv', true, 'dir', fullfile(spe11_folder, 'both_miscible'));
writeFluidFlowerPROPS(tab_h2o, tab_co2, tab_sol, 'rs', true, 'rv', false, 'dir', fullfile(spe11_folder, 'co2_miscible'));
%%
writeFluidFlowerPROPS(tab_h2o, tab_co2, tab_sol, 'rs', false, 'rv', false, 'dir', fullfile(spe11_folder, 'immiscible'));
%%
% Now the rel. perm.
writeFluidFlowerSGOF('dir', spe11_folder)
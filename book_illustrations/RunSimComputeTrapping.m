function [volume, leaked_boundary, structural_utilized, ...
            states, categorized_vols] = RunSimComputeTrapping(grid, rock, rate, state, model, ...
                                                                trapped_imperm, trapped_lowperm, swr, sor, ...
                                                                dt, bc, inj_stop_ratio)
% Run simulations for different injection rates, and return max injection
% rate that does not yield any leakage.
% Return total volume that can be injected without leakage, and what open
% boundary (top or right) caused the leakage
        
n_steps = numel(dt);
x = grid.X;
z = grid.Z;
right_cells = grid.G.cells.indexMap(x == max(x));
top_cells = grid.G.cells.indexMap(z == min(z));

dims = grid.G.cartDims;
nz = dims(3);
cc_cf_x = max(diff(x))/2; % distance from cell centroids to face centroids in x-direction
cc_cf_z = max(diff(z))/2; % NB: assumes uniform spacing!
lx = max(x) + cc_cf_x;
lz = max(z) + cc_cf_z;

leaked_boundary = 'none';

categorized_vols = {{}, {}, {}, {}}; % residual, structural, free plume, leaked

t = 0; 

well_h = 1; % cell perforations in vertical direction
perforation_idx = grid.G.cells.indexMap(z < max(z) & z >= max(z)-well_h*lz/nz & x < lx/5);
W = addWell([], grid.G, rock, perforation_idx, ...
        'Type', 'rate', 'Val', rate, ...
        'InnerProduct', 'ip_tpf', ...
        'Radius', 0.1, 'Dir', 'x', ...
        'Comp_i', [0, 1], 'Sign', 1, ... % inject CO2
        'Name', 'P1');  

schedule = simpleSchedule(dt, 'W', W, 'bc', bc);

times = cumsum(dt)/year();
%inj_stop = fix(inj_stop_ratio*n_steps);
[~, inj_stop] = min(abs(times - inj_stop_ratio*times(end)));
schedule.control(2) = schedule.control(1); % create second well
schedule.control(2).W.status = 0; % shut off second well
schedule.step.control(inj_stop:n_steps) = 2; % swap active well from first to second at halfway      

[wellSols, states] = simulateScheduleAD(state, model, schedule, 'Verbose', false);

tot_vols = zeros(numel(states)+1, 1);% first element is initial state 
simulation_vols = zeros(size(tot_vols));

structural_vols = struct;
structural_vols.imperm = zeros(size(tot_vols)); % trapped under imperm layers (excluding residual part)
structural_vols.lowperm = zeros(size(tot_vols)); % trapped under lowperm layers (excluding residual part)
structural_vols.res = zeros(size(tot_vols)); % residual part structurally trapped in IMPERM layers

residual_vols = zeros(size(tot_vols));    
free_vols = zeros(size(tot_vols));
leaked_vols = zeros(size(tot_vols));   
leaked_ratios = zeros(size(tot_vols));

for i=1:numel(states)       
    t = t + dt(i);
    Sw = states{i}.s(:,1);
    Sn = 1 - Sw;
    assert(max(Sw) < 1+eps && min(Sw) > -eps);               

    tot_vols(i+1) = sum(grid.G.cells.volumes.*Sn)*mean(rock.poro);
    simulation_vols(i+1) = min(rate*t, rate*sum(dt(1:inj_stop-1)));

    if i >= inj_stop-1
        residual_vols(i+1) = VolumeTrapping.Co2ResidualTrapped(grid.G, Sn, sor, rock);
    end
    [structural_vols.imperm(i+1), structural_vols.lowperm(i+1), ...
          structural_vols.res(i+1), ~] = VolumeTrapping.Co2StructuralTrapped(grid.G, Sn, sor, trapped_imperm, trapped_lowperm, rock);        
    leaked_vols(i+1) = simulation_vols(i+1) - tot_vols(i+1);
    leaked_ratios(i+1) = leaked_vols(i+1) / simulation_vols(i+1);
end    

% linear interpolation of residual trapping (assume gradual increase in
% imbibition over time until injection stop
dt_samples = [1, inj_stop];
res_samples = residual_vols(dt_samples);
dt_interp = 1:inj_stop;
residual_vols(1:inj_stop) = interp1(dt_samples, res_samples, dt_interp, 'linear');

for i=1:numel(states) % Need to caclulate free plume AFTER interpolation of residual
    free_vols(i+1) = max(0, tot_vols(i+1) - (residual_vols(i+1) + structural_vols.imperm(i+1) + structural_vols.lowperm(i+1))); % interpolation of residual vol could cause negative value at start --> prevent this!
end

volume = simulation_vols(end);
categorized_vols{1} = residual_vols;
categorized_vols{2} = structural_vols;
categorized_vols{3} = free_vols;
categorized_vols{4} = leaked_vols;

% Calculate leakage and compute new rates
leakage = leaked_ratios(end) > 1e-6;  % avoid round-off errors / negligable "particle" leakage  

% Compute utilized percentage of (permanently) structural traps
structural_utilized = (structural_vols.imperm(end) + structural_vols.res(end)) / (sum(grid.G.cells.volumes(trapped_imperm))*mean(rock.poro)*(1-swr)) * 100;
disp('Percentage of structural traps utilized')
disp(structural_utilized)

if leakage
    if any(Sn(right_cells)) && any(Sn(top_cells))
        leaked_boundary = ['right' 'top'];
    elseif any(Sn(right_cells))
        leaked_boundary = ['right'];          
    elseif any(Sn(top_cells))
        leaked_boundary = ['top'];    
    end        
end

end


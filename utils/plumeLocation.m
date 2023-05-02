function [tip_depth, avg_depth, z_well] = plumeLocation(model, states, schedule, snr)
%Compute location of tip and average plume over time.
% PARAMETERS:
%   model    - either fine, VE or hybrid model
%   states   - computed states from model 
%               (NB1: excluding initial state!)
%               (NB2: reconstructed fine states for VE/hybrid!)
%   schedule - schedule used for simulation
%   snr      - residual CO2 saturation
% RETURNS:
%   tip_speed - migration speed of tip of mobile plume for each state
%   avg_speed - migration speed of average position of mobile plume
%   z_well    - mean vertical location of injection well(s)

if isfield(model.G, 'parent')
    G = model.G.parent;
else
    G = model.G;   
end
[ii, jj, kk] = gridLogicalIndices(G);

num_states = numel(states);
W = schedule.control(1).W;
t = schedule.step.val;
inj_rate = schedule.control(1).W.val;
well_on = schedule.step.control == 1;
cum_inj = cumsum(t.*inj_rate.*well_on); % DOES NOT WORK?
t = cumsum(schedule.step.val);

z_well = mean(G.cells.centroids(W.cells,3)); % assume CO2 is injected uniformly through cell perforation    
k_well = fix(mean(kk(W.cells)));
k_min = min(kk);
z_min = min(G.faces.centroids(:,3));

tip_depth = zeros(num_states+1,1);
avg_depth = zeros(num_states+1,1);
tip_depth(1) = z_well;
avg_depth(1) = z_well;

reached_top = false;

for i=1:num_states
    sn = states{i}.s(:,2);    

    plume_mask = sn > snr; % 1e-5 to avoid round-off errors ?? 
    z_plume = G.cells.centroids(plume_mask, 3);

    if isempty(z_plume)
        z_tip = z_well;
        k_tip = k_well;

        z_avg = z_well;
        k_avg = k_well;
    else
        z_tip = min(z_plume);
        k_tip = min(kk(plume_mask));

        sn_scaled = sn(plume_mask) - snr; % scale by residual content so that cells with residual CO2 does no contribute to migration speed
        sn_weight = sn_scaled./sum(sn_scaled);
        z_avg = sum(z_plume.*sn_weight);
        [~, closest_idx] = min(abs(G.cells.centroids(:,3) - z_avg));
        k_avg = kk(closest_idx);
    end
    
    if (~reached_top && k_tip <= k_min) || reached_top % tip plume has reached open top boundary -> subsequent speeds unknown!
        reached_top = true;
        tip_depth(i+1) = z_min;     
    elseif t(i) > 0.01*year
        tip_depth(i+1) = min(z_tip, tip_depth(i)); % tip depth can't suddenly increase
    else
        tip_depth(i+1) = z_tip;
    end            
    
    if k_avg <= k_min
        avg_depth(i+1) = z_min;
    else
        %avg_depth(i+1) = abs(z_avg - z_well)/t;   
        avg_depth(i+1) = z_avg;
    end
end
end


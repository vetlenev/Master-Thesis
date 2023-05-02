function [tip_speed, avg_speed] = migrationSpeed(model, states, schedule, snr)
%MIGRATIONSPEED Summary of this function goes here
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

if isfield(model.G, 'parent')
    G = model.G.parent;
else
    G = model.G;   
end
[ii, jj, kk] = gridLogicalIndices(G);

t = 0;
num_states = numel(states);
W = schedule.control(1).W;
z_well = mean(G.cells.centroids(W.cells,3)); % assume CO2 is injected uniformly through cell perforation    
k_well = fix(mean(kk(W.cells)));
k_min = min(kk);

tip_speed = zeros(num_states+1,1);
avg_speed = zeros(num_states+1,1);
tip_speed(1) = 0;
avg_speed(1) = 0;

for i=1:num_states
    t = t + (schedule.step.val(i))/year();
    sn = states{i}.s(:,2);    
    
%     ii_nz = ii(sn > snr); jj_nz = jj(sn > snr);
%     kk_nz = kk(sn > snr);
%     mask_plume = intersect(find(ii_nz), ...
%                            intersect(find(jj_nz), find(kk_nz)) ...
%                            );
    %z_plume = G.cells.centroids(mask_plume, 3); 
    plume_mask = sn > 1e-5; % to avoid round-off errors % sn > snr doesn't work if open top boundary since hybrid model will never reach sn > snr for topmost VE layer
    z_plume = G.cells.centroids(plume_mask, 3);
      
    if isempty(z_plume)
        z_tip = z_well;
        k_tip = k_well;
    else
        z_tip = min(z_plume);
        k_tip = min(kk(plume_mask));
    end    
    
    if k_tip <= k_min % tip plume has reached open top boundary -> subsequent speeds unknown!
        tip_speed(i+1) = nan;
    else
        tip_speed(i+1) = abs(z_tip - z_well)/t;
    end    
        
    if i == 200
        test = 0;
    end
 
    if isempty(z_plume)
        z_avg = z_well;
        k_avg = k_well;
    else
        sn_scaled = sn(plume_mask) - snr; % scale by residual content so that cells with residual CO2 does no contribute to migration speed
        sn_weight = sn_scaled./sum(sn_scaled);
        z_avg = sum(z_plume.*sn_weight);
        [~, closest_idx] = min(abs(G.cells.centroids(:,3) - z_avg));
        k_avg = kk(closest_idx);
    end

    % Does it make sense to check this when snr is included in average
    % plume?
%     if k_avg <= k_min % avg plume has reached open top boundary -> subsequent speeds unknown!
%         avg_speed(i+1) = nan;
%     else
%         avg_speed(i+1) = abs(z_avg - z_well)/t;
%     end    
    avg_speed(i+1) = abs(z_avg - z_well)/t;
end
end


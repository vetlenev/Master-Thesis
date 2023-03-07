function schedule = getSchedule_bo(timesteps, W, bc, t)

schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', bc);                 % Simple schedule, this would inject for the total simTime
tmp = cell(2, 1);                                                           % create 2 schedules
 
schedule = struct('step', schedule_inj.step);                               % timesteps and wells for each timestep
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);                 % add 2 fields for 2 wells
schedule.control(1).W = W;                                                  % field 1 used during injection
schedule.control(2).W = W;                                                  % nr of wells must be the same for the entire simulation
schedule.control(2).W.val = 0;                                              % field 2 rate 0 (after injection)

schedule.control(1).bc = bc;
schedule.control(2).bc = bc;

schedule.step.control(cumsum(schedule.step.val) > t(1)) = 2;                % set timesteps after injTime to use Well field 2
    
end
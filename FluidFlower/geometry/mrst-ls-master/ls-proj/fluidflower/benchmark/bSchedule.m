function schedule = bSchedule(timesteps, W, bc, inj_case, opt)
%
%
%
ts = timesteps;
tsm = (cumsum(ts)/minute)';
schedule_inj = simpleSchedule(ts, 'W', W, 'bc', bc);  
%schedule = schedule_inj;

if inj_case == 1 || inj_case == 2
    injrates = [opt.rate{:, 1}]*(milli*litre)/minute;
    trd1 = opt.schedule(2)*minute + [opt.rate{2,2}];
    tru2 = opt.schedule(3)*minute + [opt.rate{3,2}];
    trd2 = opt.schedule(4)*minute + [opt.rate{4,2}];  
    tru3 = opt.schedule(5)*minute + [opt.rate{5,2}];
    trd3 = opt.schedule(6)*minute + [opt.rate{6,2}];  
    tinjr = [opt.rate{1,2} trd1 tru2 trd2 tru3 trd3]/minute;
    
    tmp = cell(numel(injrates), 1);
    schedule = struct('step', schedule_inj.step);
    schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
    for n=1:numel(injrates)
        schedule.control(n).W = W;
        schedule.control(n).bc = bc;
        schedule.control(n).W(1).val = injrates(n);
        idti = find(tsm == tinjr(n));
        if n < numel(injrates)
            idtf = find(tsm == tinjr(n+1));
        else
            idtf = numel(tsm);
        end
        schedule.step.control(idti:idtf) = n;
    end
    
elseif inj_case == 3
    injrates = [opt.rate{1} opt.rate{3,1} opt.rate{2,1}];
    tru2 = opt.schedule(3)*minute + [opt.rate{3,2}];
    if strcmp(opt.rampdown, 'fast')
        trd = opt.schedule(4)*minute + [opt.rate{4,2}];  
    else
        trd = (opt.schedule(4) + [opt.rate{4,2}])*minute; 
    end
    tinjr = [opt.rate{1,2} tru2 trd]/minute;
    
    tmp = cell(numel(injrates), 1);
    schedule = struct('step', schedule_inj.step);
    schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
    id2 = find(diff(injrates) < 0, 1) + 1;
    idrd = find(diff(injrates) < 0);    
    idrd = idrd(2)+1;
    for n=1:numel(injrates)
        schedule.control(n).W = W;
        schedule.control(n).bc = bc;
        if n < id2
            schedule.control(n).W(1).val = injrates(n);
            schedule.control(n).W(2).val = 0;
            schedule.control(n).W(2).status = false;
        elseif n < idrd
            schedule.control(n).W(1).val = W(1).val;            % max
            schedule.control(n).W(2).val = injrates(n);
            schedule.control(n).W(2).status = true;
        else
            schedule.control(n).W(1).val = injrates(n);
            schedule.control(n).W(2).val = injrates(n);
            if n == numel(injrates)
                schedule.control(n).W(1).status = false;
                schedule.control(n).W(2).status = false;
            end
            %schedule.control(n).W(2).val = 0;  % no injection in 2nd well
        end
        idti = find(tsm == tinjr(n));
        if n < numel(injrates)
            idtf = find(tsm == tinjr(n+1));
        else
            idtf = numel(tsm);
        end
        schedule.step.control(idti:idtf) = n;
    end
end

end
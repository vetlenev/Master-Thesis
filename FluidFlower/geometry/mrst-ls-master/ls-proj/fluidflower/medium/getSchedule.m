function schedule = getSchedule(inj_type, timesteps, W, bc, opt)
%
%
%
ts = timesteps;
tsm = [0 cumsum(ts)/minute]';
schedule_inj = simpleSchedule(ts, 'W', W, 'bc', bc);                 % this would inject always
injrates = [opt.rate{:, 1}]*(milli*litre)/minute;
idz = find([opt.rate{:,1}]==0);
tmp = cell(numel(injrates), 1);
schedule = struct('step', schedule_inj.step);
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
tru1 = [0 opt.rate{1,2}];
if inj_type == 1 
    tru2 = [opt.rate{3,2} 6];
    trd2 = opt.rate{end,2}/minute;
    trd1 = fliplr(opt.schedule{2}-[0 opt.rate{2,2}]);
    tru2 = opt.schedule{3}(1)+(opt.schedule{3}(2)/minute) + ...
                     (tru2-1);
    trd2 = fliplr(opt.schedule{4}-cumsum([0 fliplr(trd2)]));             
elseif inj_type == 2
    trd1 = opt.rate{2,2};
    trd1 = (fliplr(opt.schedule{2}(1)*minute + opt.schedule{2}(2) ...
           -cumsum([0 fliplr(trd1)])))/minute; 
    tru2 = [];
    trd2 = [];
elseif inj_type == 3
    trd1 = fliplr(opt.schedule{2}(1)-[0 opt.rate{2,2}] + opt.schedule{2}(2)/minute);
    tru2 = [opt.rate{3,2} 4];
    tru2 = opt.schedule{3}(1)+(opt.schedule{3}(2)/minute) + (tru2-1);
    trd2 = fliplr(opt.schedule{4}(1)-[0 opt.rate{4,2}] + opt.schedule{4}(2)/minute);
end
tinjr = [tru1 trd1 tru2 trd2];

for n=1:numel(injrates)
    schedule.control(n).W = W; 
    schedule.control(n).bc = bc;
    if inj_type ~= 2
        if n <= idz(1)
            schedule.control(n).W(1).val = injrates(n);
            schedule.control(n).W(2).val = 0;
        else
            schedule.control(n).W(1).val = 0;
            schedule.control(n).W(2).val = injrates(n);
            %schedule.control(n).W(2).val = 0;  % no injection in 2nd well
        end
    else
         schedule.control(n).W.val = injrates(n);
    end
    idti = find(tsm == tinjr(n));
    if n < numel(injrates)
        idtf = find(tsm == tinjr(n+1));
    else
        idtf = numel(tsm)-1;
    end
    schedule.step.control(idti:idtf) = n;
end 

end
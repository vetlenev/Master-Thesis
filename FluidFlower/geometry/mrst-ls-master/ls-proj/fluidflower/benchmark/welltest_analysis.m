% Load simulation and experiment data, plot results, pick best parameter
% combination based on MAE.

clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils 
mrstVerbose off

% Simulation opts
load('base_problem_data2_wt2.mat');  % this just for loading output, if rock,etc
                                     % needs to be viewed we need to compute 
                                     % rock,etc for each case.
topDir = 'D:\Lluis\fluidflower\benchmark\welltest\wt2';
%batch = 'all2'; 
%[~, ~, m] = bSimOpts(inj_case, mesh_size, [], false, batch);
nsim = size(m, 1);
%nsim = 100;

fports = [0.0120 0.53 1.01; 0.0120 0.53 0.61];
dist = pdist2(G.cells.centroids, fports);
[~, fpInx] = min(dist,[],1);
wellId = [wellId; fpInx'];

% Simulation data
%wellSols = cell(nsim, 1);
%states = cell(nsim, 1);
dp = cell(nsim, 1);
bhp = zeros(numel(timesteps), nsim);
wellId_dp = cell(nsim, 1);
wellId_dbhp = zeros(numel(timesteps), nsim);
tic
parfor k=1:nsim           % careful with memory!
    folderName = ['exp' num2str(inj_case) '_bmesh' num2str(mesh_size) ...
                  '_modelcase' num2str(model_case) ...
                  '_reportedProps_stepMaxRate_kmult_ESF', num2str(m(k,1)) ...
                  '_ESFsup', num2str(m(k,2)), '_C', num2str(m(k,3)), '_D', num2str(m(k,4)) ...
                  '_E', num2str(m(k,5)), '_F', num2str(m(k,6)), '_G', num2str(m(k,7))];
              
    problem = packSimulationProblem(state0, model, schedule, folderName, ...
                                    'Name', folderName, 'Directory', ...
                                     topDir, 'NonLinearSolver', nls);
    [wellSols, states] = getPackedSimulatorOutput(problem);
    
    dp{k} = reshape(cell2mat(cellfun(@(x) x.pressure - state0.pressure, ...
                    states, 'UniformOutput', false)), G.cells.num, numel(timesteps));
    bhp(:, k) = cellfun(@(x) x.bhp, wellSols);
    wellId_dbhp(:, k) = cellfun(@(x) x.bhp - state0.pressure(W.cells), wellSols);
    wellId_dp{k} = reshape(cell2mat(cellfun(@(x) x.pressure(wellId) - state0.pressure(wellId), ...
                                    states, 'UniformOutput', false)), ...
                           numel(wellId), numel(timesteps))';
    
    if mod(k, 50) == 0
        disp(['Loaded ' num2str(k) ' / ' num2str(nsim)])
    end
    
end
toc

% Experimental data
pth = 'C:\Users\Lluis\matlab\mrst-dev\mrst-ls\ls-proj\fluidflower\benchmark\welltests_data\';
np = 9;
if inj_case == 1
    fn = 'pressure_measurements_welltest_port_17_7.xls';
    id_tcol = [7 7+(3:3:3*(np-1))];
    id_pcol = [8 8+(3:3:3*(np-1))];
    id_t0 = [repelem(1143, 6) 1142, 1143, 1142];
elseif inj_case == 2
    fn = 'pressure_measurements_welltest_port_9_3.xls';
    id_tcol = [4 4+(3:3:3*np)]; id_tcol(id_tcol==28) = [];
    id_pcol = [5 5+(3:3:3*np)]; id_pcol(id_pcol==29) = [];
    id_t0 = repelem(228, 9);
end
data = readtable([pth, fn], 'Sheet', 1);
p_data = cell(1, np);
t_data = cell(1, np);
ts = cumsum(timesteps);
for n=1:np
    tval = table2array(data(:, id_tcol(n)));
    t0 = tval(id_t0(n));
    tval = tval - t0;
    tfin = opt.schedule(end)*minute;
    [~, id_tfin] = min(abs(tfin - tval));
    t_data{n} = tval(id_t0(n):id_tfin);
    p_data{n} = table2array(data(id_t0(n):id_tfin, id_pcol(n)))*1e3;    % mbar
    [tip, tid] = unique(t_data{n});
    for j = 1:numel(timesteps)
        [~, id_t] = min(abs(t_data{n} - ts(j)));
        i1 = max([1 id_t-5]);
        i2 = min([numel(t_data{n}) id_t+5]);
        p_data_ip(j,n) = mean(p_data{n}(i1:i2));
    end
%    p_data_ip(:, n) = interp1(tip, p_data{n}(tid), cumsum(timesteps), 'linear');
end
if inj_case == 1
    id_i1 = 4;
    id_i2 = 1;
    id_p1 = [5 7];
    id_p2 = [8 9];
    id_f1 = [2 6];
    id_f2 = 3;
elseif inj_case == 2
    id_i1 = 4;
    id_i2 = 8;
    id_p1 = [1 9];
    id_p2 = [5 7];
    id_f1 = [2 6];
    id_f2 = 3;
end
p_i2 = p_data{id_i2};            % 0-2.5 bar
p_i1 = p_data{id_i1};            % -1-2.5 bar
p_p1_1 = p_data{id_p1(1)};       % 0-6 bar
p_p1_2 = p_data{id_p1(2)};       % 0-4 bar
p_p2_1 = p_data{id_p2(1)};       % 0-100,
p_p2_2 = p_data{id_p2(2)};       % 0-4 bar
p_f1_1 = p_data{id_f1(1)};       % 0-2.5,
p_f1_2 = p_data{id_f1(2)};       % 0-100 bar
p_f2 = p_data{id_f2};            % 0-2.5 bar

p_i2_ip = p_data_ip(:,id_i2);           
p_i1_ip = p_data_ip(:,id_i1);          
p_p1_1_ip = p_data_ip(:,id_p1(1));       
p_p1_2_ip = p_data_ip(:,id_p1(2));       
p_p2_1_ip = p_data_ip(:,id_p2(1));       
p_p2_2_ip = p_data_ip(:,id_p2(2));       
p_f1_1_ip = p_data_ip(:,id_f1(1));       
p_f1_2_ip = p_data_ip(:,id_f1(2));       
p_f2_ip = p_data_ip(:,id_f2);            

% Separate
wellId_dp_1 = reshape(cell2mat(cellfun(@(x) x(:,1), wellId_dp, 'UniformOutput', false)), numel(timesteps), nsim);
wellId_dp_2 = reshape(cell2mat(cellfun(@(x) x(:,2), wellId_dp, 'UniformOutput', false)), numel(timesteps), nsim);
wellId_dp_3 = reshape(cell2mat(cellfun(@(x) x(:,3), wellId_dp, 'UniformOutput', false)), numel(timesteps), nsim);
wellId_dp_4 = reshape(cell2mat(cellfun(@(x) x(:,4), wellId_dp, 'UniformOutput', false)), numel(timesteps), nsim);
wellId_dp_5 = reshape(cell2mat(cellfun(@(x) x(:,5), wellId_dp, 'UniformOutput', false)), numel(timesteps), nsim);
wellId_dp_6 = reshape(cell2mat(cellfun(@(x) x(:,6), wellId_dp, 'UniformOutput', false)), numel(timesteps), nsim);

% Plot
latx={'Interpreter','latex'};
prpl = [92, 16, 199]/255;
figure(37)
subplot(2,3,1)
    hold on
    p1 = plot(cumsum(timesteps)/minute, wellId_dp_1/100, '-', ...
              'color', [0.8 0.8 0.8], 'linewidth', 1, ...
              'DisplayName', 'Simulation');
    p2 = plot(t_data{id_i1}/minute, p_i1, '-k', 'linewidth', 1, ...
              'DisplayName', '[-1--2.5] bar');
    p3 = plot(cumsum(timesteps)/minute, p_i1_ip, '-b', 'linewidth', 0.5, ...
              'DisplayName', 'Interpolated');
         
%     if inj_case == 1
  
%     elseif inj_case == 2
%         p1 = plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dbhp/100, '-r', 'linewidth', 1, ...
%                      'DisplayName', 'Simulation');
%     end
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('I$_1$ $\vert$ $y, z = 0.93, 1.01$m $\vert$ [9, 3]', latx{:})
    grid on
    xlabel('$t$ [min]', latx{:})
    ylabel('$\Delta p$ [mbar]', latx{:})
    h=legend([p1(1), p2, p3], 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,2)
    hold on
    plot(cumsum(timesteps)/minute, wellId_dp_4/100, '-', ...
            'color', [0.8 0.8 0.8], 'linewidth', 1)
    p12 = plot(t_data{id_i2}/minute, p_i2, '-k', 'linewidth', 1, ...
         'DisplayName', '[0--2.5] bar');
    plot(cumsum(timesteps)/minute, p_i2_ip, '-b', 'linewidth', 0.5, ...
         'DisplayName', 'Interpolated');
%     if inj_case == 2
        
%     elseif inj_case == 1
%         plot(cumsum(ts_out{1}.timesteps)/minute, wellId_dbhp/100, '-r', 'linewidth', 1)
%     end
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('I$_2$ $\vert$ $y, z = 1.73, 0.61$m $\vert$ [17, 7]', latx{:})
    grid on
    
subplot(2,3,3)
    hold on
    plot(cumsum(timesteps)/minute, wellId_dp_2/100, '-', ...
        'color', [0.8 0.8 0.8], 'linewidth', 1)
    p13 = plot(t_data{id_p1(1)}/minute, p_p1_1, '-k', 'linewidth', 1, ...
         'DisplayName', '[0--6] bar');
    p23 = plot(t_data{id_p1(2)}/minute, p_p1_2, '-', 'color', [0.4 0.4 0.4], ...
         'linewidth', 1, 'DisplayName', '[0--4] bar');
    plot(cumsum(timesteps)/minute, p_p1_1_ip, '-b', 'linewidth', 0.5, ...
         'DisplayName', 'Interpolated');
    plot(cumsum(timesteps)/minute, p_p1_2_ip, '-c', 'linewidth', 0.5, ...
         'DisplayName', 'Interpolated');
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('P$_1$ $\vert$ $y, z = 1.53, 0.81$m $\vert$ [15, 5]', latx{:})
    grid on
    
subplot(2,3,4)
    hold on
    plot(cumsum(timesteps)/minute, wellId_dp_3/100, '-', ...
         'color', [0.8 0.8 0.8], 'linewidth', 1)
    p14 = plot(t_data{id_p2(1)}/minute, p_p2_1, '-k', 'linewidth', 1, ...
         'DisplayName', '[0--100] bar');
    p24 = plot(t_data{id_p2(2)}/minute, p_p2_2, '-', 'color', [0.4 0.4 0.4], ...
         'linewidth', 1, 'DisplayName', '[0--4] bar');
    plot(cumsum(timesteps)/minute, p_p2_1_ip, '-b', 'linewidth', 0.5, ...
         'DisplayName', 'Interpolated');
    plot(cumsum(timesteps)/minute, p_p2_2_ip, '-c', 'linewidth', 0.5, ...
         'DisplayName', 'Interpolated');
    grid on
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('P$_2$ $\vert$ $y, z = 1.73, 0.21$m $\vert$ [17, 11]', latx{:})
    hold off

subplot(2,3,5)
    plot(cumsum(timesteps)/minute, wellId_dp_5/100, '-', ...
         'color', [0.8 0.8 0.8], 'linewidth', 1)
    hold on
    p15 = plot(t_data{id_f1(1)}/minute, p_f1_1, '-k', 'linewidth', 1, ...
         'DisplayName', '[0--2.5] bar');
    p25 = plot(t_data{id_f1(2)}/minute, p_f1_2, '-', 'color', [0.6 0.6 0.6], ...
         'linewidth', 1, 'DisplayName', '[0--100] bar');
    plot(cumsum(timesteps)/minute, p_f1_1_ip, '-b', 'linewidth', 0.5, ...
         'DisplayName', 'Interpolated');
    plot(cumsum(timesteps)/minute, p_f1_2_ip, '-c', 'linewidth', 0.5, ...
         'DisplayName', 'Interpolated');
    grid on
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('F$_1$ $\vert$ $y, z = 0.53, 1.01$m $\vert$ [5, 3]', latx{:})
    hold off
    
subplot(2,3,6)
    hold on
    plot(cumsum(timesteps)/minute, wellId_dp_6/100, '-', ...
         'color', [0.8 0.8 0.8], 'linewidth', 1)
    p16 = plot(t_data{id_f2}/minute, p_f2, '-k', 'linewidth', 1, ...
         'DisplayName', '[0--2.5] bar');
    plot(cumsum(timesteps)/minute, p_f2_ip, '-b', 'linewidth', 0.5, ...
         'DisplayName', 'Interpolated');
    grid on
    %ylim([29 31])
    xlim([0 180])
    xticks([0 30 60 90 120 150 180])
    title('F$_2$ $\vert$ $y, z = 0.53, 0.61$m $\vert$ [5, 7]', latx{:})
    hold off

% Best Parameter combination
mae = zeros(nsim, np-4);
for n=1:nsim
    mae(n, 1) = sum(abs(wellId_dp_1(:,n)/100 - p_i1_ip))/numel(ts);
    %mae(n, 2) = sum(abs(wellId_dp_4(:,n)/100 - p_i2_ip))/numel(ts);
    mae(n, 2) = sum(abs(wellId_dp_2(:,n)/100 - p_p1_2_ip))/numel(ts);
    mae(n, 3) = sum(abs(wellId_dp_3(:,n)/100 - p_p2_2_ip))/numel(ts);
    mae(n, 4) = sum(abs(wellId_dp_5(:,n)/100 - p_f1_1_ip))/numel(ts);
    mae(n, 5) = sum(abs(wellId_dp_6(:,n)/100 - p_f2_ip))/numel(ts);
end
mae_sum = sum(mae, 2);
mae_std = std(mae, [], 2);
[~, idmin] = mink(mae_sum, 1);
std_idmin = mae_std(idmin);
[~, idmin_std] = mink(mae_std, 1);
mae_idmin_std = mae_sum(idmin_std);


% Plot 
figure(99)
errorbar(1:nsim, mae_sum, mae_std, 's', 'color', 'k')
hold on
errorbar(idmin, mae_sum(idmin), std_idmin, 's', 'color', 'g')

figure(37)
subplot(2,3,1)
hold on
plot(cumsum(timesteps)/minute, wellId_dp_1(:,idmin)/100, '-g', ...
     'linewidth', 1, 'DisplayName', 'Min MAE');
plot(cumsum(timesteps)/minute, wellId_dp_1(:,idmin_std)/100, '-m', ...
     'linewidth', 1, 'DisplayName', 'Min std');
subplot(2,3,2)
hold on
plot(cumsum(timesteps)/minute, wellId_dp_4(:,idmin)/100, '-g','linewidth', 1)
plot(cumsum(timesteps)/minute, wellId_dp_4(:,idmin_std)/100, '-m','linewidth', 1)
h=legend(p12, 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,3)
hold on
plot(cumsum(timesteps)/minute, wellId_dp_2(:,idmin)/100, '-g', 'linewidth', 1)
plot(cumsum(timesteps)/minute, wellId_dp_2(:,idmin_std)/100, '-m','linewidth', 1)
h=legend([p13, p23], 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,4)
hold on
plot(cumsum(timesteps)/minute, wellId_dp_3(:,idmin)/100, '-g', 'linewidth', 1)
plot(cumsum(timesteps)/minute, wellId_dp_3(:,idmin_std)/100, '-m','linewidth', 1)
h=legend([p14, p24], 'location', 'northeast', latx{:}, 'fontsize', 8);
set(h.BoxFace, 'ColorType','truecoloralpha', ...
    'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,5)
hold on
plot(cumsum(timesteps)/minute, wellId_dp_5(:,idmin)/100, '-g', 'linewidth', 1)
plot(cumsum(timesteps)/minute, wellId_dp_5(:,idmin_std)/100, '-m','linewidth', 1)
h=legend([p15, p25], 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
                'ColorData', uint8(255*[1;1;1;.5]));
subplot(2,3,6)
hold on
plot(cumsum(timesteps)/minute, wellId_dp_6(:,idmin)/100, '-g', 'linewidth', 1)
plot(cumsum(timesteps)/minute, wellId_dp_6(:,idmin_std)/100, '-m','linewidth', 1)
h=legend(p16, 'location', 'northeast', latx{:}, 'fontsize', 8);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
        'ColorData', uint8(255*[1;1;1;.5]));
    
    
% Val min
% ESF, ESFsup, C, D, E, F, G
k_match(1,:) = m(idmin,:);
k_match(2,:) = m(idmin_std,:);

%% Plot smear-based fault permeability algorithm convergence results
% Save data individually from faultPermSmearThrow_base. Here we load and
% plot perm as a function of cell resolution (length) in the vertical
% direction, since it's the one controlling the square matrix and grid
% dimensions.

clc, clear, close all force 
mrstModule add ls-proj ls-utils coarsegrid upscaling incomp mimetic mpfa ...
           streamlines mrst-gui

% Load data for different stratigraphies
dir = 'C:\Users\lsalo\matlab\sim_data\fault_perm_algorithm\singleParamVariations\nSim_analysis\';
fname = {'nSimSmear_2beds_fT2fdip60_phi1030_SSFc1.5_LsSmear', ...
         'nSimSmear_2beds_fT2fdip60_phi1030_SSFc1.8_LsSmear', ...
         'nSimSmear_2beds_fT4fdip60_phi1030_SSFc1.5_Ls50', ...
         'nSimSmear_2beds_fT4fdip60_phi1030_SSFc1.8_Ls50', ...
         'nSimSmear_4beds_fT2fdip60_phi1030_SSFc2.5_LsSmear', ...
         'nSimSmear_4beds_fT2fdip60_phi1030_SSFc3_LsSmear', ...
         'nSimSmear_4beds_fT2fdip60_phi1030_SSFc3.5_LsSmear', ...
         'nSimSmear_6beds_fT1fdip60_phi1030_SSFc2_Ls20', ...
         'nSimSmear_6beds_fT1fdip60_phi1030_SSFc5_LsSmear', ...
         'nSimSmear_8beds_fT2fdip60_phi530_SSFc7_Ls30', ...
         'nSimSmear_10beds_fT1fdip60_phi1030_SSFc7_Ls40', ...
         'nSimSmear_10beds_fT5fdip60_phi530_SSFc6_Ls10', ...
         'nSimSmear_20beds_fT0.5fdip60_phi535_SSFc10_LsSmear', ...
         'nSimSmear_20beds_fT1fdip60_phi535_SSFc10_LsSmear', ...
         'nSimSmear_6_8bedsIrreg_fT1fdip60_phi1030_SSFc3_LsSmear', ...
         'nSimSmear_6_8bedsIrreg_fT2fdip60_phi1030_SSFc5_LsSmear'};

for n=1:numel(fname)
   data = load([dir fname{n} '.mat']);
   
   perms = cell2mat(data.permMed);
   permx = perms(:, 1);
   permy = perms(:, 2);
   permz = perms(:, 3);
   pctErr = cell2mat(data.permPctErr);
   pct10x = permx - pctErr(:,1);
   pct10y = permy - pctErr(:,2);
   pct10z = permz - pctErr(:,3);
   pct90x = permx + pctErr(:,4);
   pct90y = permy + pctErr(:,5);
   pct90z = permz + pctErr(:,6);
   
   rx(:, n) = abs(log10(permx(1:end-1)) - log10(permx(end)));
   ry(:, n) = abs(log10(permy(1:end-1)) - log10(permy(end)));
   rz(:, n) = abs(log10(permz(1:end-1)) - log10(permz(end)));
   rPct10x(:, n) = abs(log10(pct10x(1:end-1)) - log10(pct10x(end)));
   rPct10y(:, n) = abs(log10(pct10y(1:end-1)) - log10(pct10y(end)));
   rPct10z(:, n) = abs(log10(pct10z(1:end-1)) - log10(pct10z(end)));
   rPct90x(:, n) = abs(log10(pct90x(1:end-1)) - log10(pct90x(end)));
   rPct90y(:, n) = abs(log10(pct90y(1:end-1)) - log10(pct90y(end)));
   rPct90z(:, n) = abs(log10(pct90z(1:end-1)) - log10(pct90z(end)));
   
   % Plot error (with respect to metrics for 5000 sims) vs Nsim
   if n == 1
       nSim = data.nSimSmear;
       latx = {'Interpreter','latex'};
       colrsx = cmocean('thermal', 16);
       %colrsy = cmocean('matter', 8);
       colrsz = cmocean('thermal', 16);
       %colrsz(1, :) = [0.9861 0.8415 0.6061];
       
       f9 = figure(99);
       subplot(3,2,1)
       plot([500 500], [0 3.5], '-k', 'lineWidth', 1)
       hold on
       subplot(3,2,2)
       plot([500 500], [0 3.5], '-k', 'lineWidth', 1)
       hold on
       subplot(3,2,3)
       plot([500 500], [0 3.5], '-k', 'lineWidth', 1)
       hold on
       subplot(3,2,4)
       plot([500 500], [0 3.5], '-k', 'lineWidth', 1)
       hold on
       subplot(3,2,5)
       plot([500 500], [0 3.5], '-k', 'lineWidth', 1)
       hold on
       subplot(3,2,6)
       plot([500 500], [0 3.5], '-k', 'lineWidth', 1)
       hold on
   end
    
   subplot(3,2,1)
   semilogx(nSim(1:end-1), rx(:,n), '.-', 'color', colrsx(n,:))
   subplot(3,2,2)
   semilogx(nSim(1:end-1), rz(:,n), '.-', 'color', colrsz(n,:))
   subplot(3,2,3)
   loglog(nSim(1:end-1), rPct10x(:,n), '.-', 'color', colrsx(n,:))
   subplot(3,2,4)
   loglog(nSim(1:end-1), rPct10z(:,n), '.-', 'color', colrsz(n,:))
   subplot(3,2,5)
   loglog(nSim(1:end-1), rPct90x(:,n), '.-', 'color', colrsx(n,:))
   subplot(3,2,6)
   loglog(nSim(1:end-1), rPct90z(:,n), '.-', 'color', colrsz(n,:))
   
   if n == numel(fname)
       subplot(3,2,1)
       hold off
       xlim([5 2000])
       ylim([0 3.5])
       yticks([0 1 2 3 3.5])
       xticks([10 100 1000])
       xticklabels({'10^1', '10^2', '10^3'})
       xlabel('$N_\mathrm{sim}$', latx{:}, 'fontSize', 12)
       ylabel('$|\log \tilde{k}_{i} - \log \tilde{k}_{5000}|$ [-]', latx{:}, 'fontSize', 10)
       title('$k_{xx}$', latx{:}, 'fontSize', 13)
       grid on
       set(gca,'XScale','log')
       
       subplot(3,2,2)
       hold off
       xlim([5 2000])
       ylim([0 1.5])
       yticks([0 0.5 1 1.5])
       xticks([10 100 1000])
       xticklabels({'10^1', '10^2', '10^3'})
       xlabel('$N_\mathrm{sim}$', latx{:}, 'fontSize', 12)
       ylabel('$|\log \tilde{k}_{i} - \log \tilde{k}_{5000}|$ [-]', latx{:}, 'fontSize', 10)
       title('$k_{zz}$', latx{:}, 'fontSize', 13)
       grid on
       set(gca,'XScale','log')
       hold off
       
       subplot(3,2,3)
       xlim([5 2000])
       ylim([0 3.5])
       yticks([0 1 2 3 3.5])
       xticks([10 100 1000])
       xticklabels({'10^1', '10^2', '10^3'})
       ylabel('$|\log p_{10}k_{i} - \log p_{10}k_{5000}|$ [-]', latx{:}, 'fontSize', 10)
       grid on
       set(gca,'XScale','log')
       
       subplot(3,2,4)
       xlim([5 2000])
       ylim([0 1.5])
       yticks([0 0.5 1 1.5])
       xticks([10 100 1000])
       xticklabels({'10^1', '10^2', '10^3'})
       ylabel('$|\log p_{10}k_{i} - \log p_{10}k_{5000}|$ [-]', latx{:}, 'fontSize', 10)
       grid on
       set(gca,'XScale','log')
       
       subplot(3,2,5)
       xlim([5 2000])
       ylim([0 3.5])
       yticks([0 1 2 3 3.5])
       xticks([10 100 1000])
       xticklabels({'10^1', '10^2', '10^3'})
       ylabel('$|\log p_{90}k_{i} - \log p_{90}k_{5000}|$ [-]', latx{:}, 'fontSize', 10)
       grid on
       set(gca,'XScale','log')
       
       subplot(3,2,6)
       xlim([5 2000])
       ylim([0 1.5])
       yticks([0 0.5 1 1.5])
       xticks([10 100 1000])
       xticklabels({'10^1', '10^2', '10^3'})
       ylabel('$|\log p_{90}k_{i} - \log p_{90}k_{5000}|$ [-]', latx{:}, 'fontSize', 10)
       grid on
       set(gca,'XScale','log')
       set(f9, 'position', [500, 200, 650, 550]);
   end
   
end

% subplot(3,2,2)
% plot([200 200], [10^-3 10^2], '-k')
% hold on
% loglog(nSim(1:end-1), rz(:,1), '.-', 'color', colrsz(1,:))
% loglog(nSim(1:end-1), rz(:,2), '.-', 'color', colrsz(2,:))
% loglog(nSim(1:end-1), rz(:,3), '.-', 'color', colrsz(3,:))
% loglog(nSim(1:end-1), rz(:,4), '.-', 'color', colrsz(4,:))
% loglog(nSim(1:end-1), rz(:,5), '.-', 'color', colrsz(5,:))
% loglog(nSim(1:end-1), rz(:,6), '.-', 'color', colrsz(6,:))
% loglog(nSim(1:end-1), rz(:,7), '.-', 'color', colrsz(7,:))
% loglog(nSim(1:end-1), rz(:,8), '.-', 'color', colrsz(8,:))
% hold off
% xlim([5 2000])
% ylim([0 3])
% yticks([0 1 2 3])
% xticks([10 100 1000])
% xticklabels({'10^1', '10^2', '10^3'})
% xlabel('$N_\mathrm{sim}$', latx{:}, 'fontSize', 12)
% ylabel('$|\log k_{i} - \log k_{5000}|$ [-]', latx{:}, 'fontSize', 12)
% title('Med($k_{zz}$)', latx{:}, 'fontSize', 13)
% grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% subplot(3,2,3)
% plot([200 200], [10^-3 10^2], '-k')
% hold on
% loglog(nSim(1:end-1), rPct10x(:,1), '.-', 'color', colrsx(1,:))
% loglog(nSim(1:end-1), rPct10x(:,2), '.-', 'color', colrsx(2,:))
% loglog(nSim(1:end-1), rPct10x(:,3), '.-', 'color', colrsx(3,:))
% loglog(nSim(1:end-1), rPct10x(:,4), '.-', 'color', colrsx(4,:))
% loglog(nSim(1:end-1), rPct10x(:,5), '.-', 'color', colrsx(5,:))
% loglog(nSim(1:end-1), rPct10x(:,6), '.-', 'color', colrsx(6,:))
% loglog(nSim(1:end-1), rPct10x(:,7), '.-', 'color', colrsx(7,:))
% loglog(nSim(1:end-1), rPct10x(:,8), '.-', 'color', colrsx(8,:))
% hold off
% xlim([5 2000])
% ylim([0 3])
% yticks([0 1 2 3])
% xticks([10 100 1000])
% xticklabels({'10^1', '10^2', '10^3'})
% title('$10^{\mathrm{th}}$ pp, $k_{xx}$', latx{:}, 'fontSize', 13)
% grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% subplot(3,2,4)
% plot([200 200], [10^-3 10^2], '-k')
% hold on
% loglog(nSim(1:end-1), rPct10z(:,1), '.-', 'color', colrsz(1,:))
% loglog(nSim(1:end-1), rPct10z(:,2), '.-', 'color', colrsz(2,:))
% loglog(nSim(1:end-1), rPct10z(:,3), '.-', 'color', colrsz(3,:))
% loglog(nSim(1:end-1), rPct10z(:,4), '.-', 'color', colrsz(4,:))
% loglog(nSim(1:end-1), rPct10z(:,5), '.-', 'color', colrsz(5,:))
% loglog(nSim(1:end-1), rPct10z(:,6), '.-', 'color', colrsz(6,:))
% loglog(nSim(1:end-1), rPct10z(:,7), '.-', 'color', colrsz(7,:))
% loglog(nSim(1:end-1), rPct10z(:,8), '.-', 'color', colrsz(8,:))
% hold off
% xlim([5 2000])
% ylim([0 3])
% yticks([0 1 2 3])
% xticks([10 100 1000])
% xticklabels({'10^1', '10^2', '10^3'})
% title('$10^{\mathrm{th}}$ pp, $k_{zz}$', latx{:}, 'fontSize', 13)
% grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% subplot(3,2,5)
% plot([200 200], [10^-3 10^2], '-k')
% hold on
% loglog(nSim(1:end-1), rPct90x(:,1), '.-', 'color', colrsx(1,:))
% loglog(nSim(1:end-1), rPct90x(:,2), '.-', 'color', colrsx(2,:))
% loglog(nSim(1:end-1), rPct90x(:,3), '.-', 'color', colrsx(3,:))
% loglog(nSim(1:end-1), rPct90x(:,4), '.-', 'color', colrsx(4,:))
% loglog(nSim(1:end-1), rPct90x(:,5), '.-', 'color', colrsx(5,:))
% loglog(nSim(1:end-1), rPct90x(:,6), '.-', 'color', colrsx(6,:))
% loglog(nSim(1:end-1), rPct90x(:,7), '.-', 'color', colrsx(7,:))
% loglog(nSim(1:end-1), rPct90x(:,8), '.-', 'color', colrsx(8,:))
% hold off
% xlim([5 2000])
% ylim([0 3])
% yticks([0 1 2 3])
% xticks([10 100 1000])
% xticklabels({'10^1', '10^2', '10^3'})
% title('$90^{\mathrm{th}}$ pp, $k_{xx}$', latx{:}, 'fontSize', 13)
% grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')
% 
% subplot(3,2,6)
% plot([200 200], [10^-3 10^2], '-k')
% hold on
% loglog(nSim(1:end-1), rPct90z(:,1), '.-', 'color', colrsz(1,:))
% loglog(nSim(1:end-1), rPct90z(:,2), '.-', 'color', colrsz(2,:))
% loglog(nSim(1:end-1), rPct90z(:,3), '.-', 'color', colrsz(3,:))
% loglog(nSim(1:end-1), rPct90z(:,4), '.-', 'color', colrsz(4,:))
% loglog(nSim(1:end-1), rPct90z(:,5), '.-', 'color', colrsz(5,:))
% loglog(nSim(1:end-1), rPct90z(:,6), '.-', 'color', colrsz(6,:))
% loglog(nSim(1:end-1), rPct90z(:,7), '.-', 'color', colrsz(7,:))
% loglog(nSim(1:end-1), rPct90z(:,8), '.-', 'color', colrsz(8,:))
% hold off
% xlim([5 2000])
% ylim([0 3])
% yticks([0 1 2 3])
% xticks([10 100 1000])
% xticklabels({'10^1', '10^2', '10^3'})
% title('$90^{\mathrm{th}}$ pp, $k_{zz}$', latx{:}, 'fontSize', 13)
% grid on
% set(gca,'YScale','log')
% set(gca,'XScale','log')

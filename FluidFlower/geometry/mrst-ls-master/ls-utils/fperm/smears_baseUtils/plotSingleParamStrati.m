% Plot single parameter variations for several strati in one plot
%
%
%% Plot smear-based fault permeability algorithm convergence results
% Save data individually from faultPermSmearThrow_base. Here we load and
% plot perm as a function of cell resolution (length) in the vertical
% direction, since it's the one controlling the square matrix and grid
% dimensions.

clc, clear, close all force 

% Load data for different stratigraphies
plotWhat = 'SSFc';
switch plotWhat
    case 'clayResidual'
dir = 'C:\Users\lsalo\matlab\sim_data\fault_perm_algorithm\singleParamVariations\clayResidual_strati\';
fname = {'clayResidual_2beds_fT1fdip60_phis30_SSFc5_LsSmear', ...
         'clayResidual_4beds_fT1fdip60_phis30_SSFc5_LsSmear', ...
         'clayResidual_5beds_fT1fdip60_phis30_SSFc5_LsSmear', ...
         'clayResidual_6beds_fT1fdip60_phis30_SSFc5_LsSmear', ...
         'clayResidual_8beds_fT1fdip60_phis30_SSFc5_LsSmear', ...
         'clayResidual_10beds_fT1fdip60_phis30_SSFc5_LsSmear'};
    case 'SSFc'
        dir = 'C:\Users\lsalo\matlab\sim_data\fault_perm_algorithm\singleParamVariations\SSFc_strati\';
        fname = {'SSFc_2beds_fT1fdip60_phi1030_LsSmear', ...
                 'SSFc_4beds_fT1fdip60_phi1030_LsSmear', ...
                 'SSFc_8beds_fT1fdip60_phi1030_LsSmear', ...
                 'SSFc_10beds_fT1fdip60_phi1030_LsSmear'};
end

permMed = cell(numel(fname), 1);
for n=1:numel(fname)
   data = load([dir fname{n} '.mat']);
   permMed{n} = reshape([data.Perm.estimate], 3, numel(data.Perm))';
   %PermPctErr = reshape([data.Perm.pctErr], 6, numel(data.Perm))';
   %pct10{n} = permMed{n} - PermPctErr(:,1:3);
   %pct90{n} = permMed{n} + PermPctErr(:,4:6);
   
   switch plotWhat
       case 'clayResidual'
           if n == 1
               latx = {'Interpreter','latex'};
               colrs = winter(numel(fname));
               
               f9 = figure(99);
               cr = data.clayResidual;
               subplot(1,3,1)
               text(14, 1.5e-3, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
               hold on
               text(14, 1.5e-2, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
               text(12, 150,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
               semilogy([cr(1) cr(end)], [data.claySandPerms(1) data.claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               semilogy([cr(1) cr(end)], [data.claySandPerms(2) data.claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               semilogy([cr(1) cr(end)], [data.claySandPerms(3) data.claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               
               subplot(1,3,2)
               %text(14, 1.5e-3, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
               hold on
               %text(14, 1.5e-2, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
               %text(12, 150,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
               semilogy([cr(1) cr(end)], [data.claySandPerms(1) data.claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               semilogy([cr(1) cr(end)], [data.claySandPerms(2) data.claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               semilogy([cr(1) cr(end)], [data.claySandPerms(3) data.claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               
               subplot(1,3,3)
               %text(14, 1.5e-3, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
               hold on
               %text(14, 1.5e-2, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
               %text(12, 150,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
               semilogy([cr(1) cr(end)], [data.claySandPerms(1) data.claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               semilogy([cr(1) cr(end)], [data.claySandPerms(2) data.claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               semilogy([cr(1) cr(end)], [data.claySandPerms(3) data.claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)              
           end
           
           subplot(1,3,1)
           p1 = semilogy(cr, permMed{n}(:,1), '.-', 'color', colrs(n, :));
           subplot(1,3,2)
           p2 = semilogy(cr, permMed{n}(:,2), '.-', 'color', colrs(n, :));
           subplot(1,3,3)
           p3 = semilogy(cr, permMed{n}(:,3), '.-', 'color', colrs(n, :));
           
           if n == numel(fname)
               subplot(1,3,1)
               hold off
               xlabel('$\phi^`_\mathrm{r, c}$', 'Interpreter', 'latex', 'fontSize', 13)
               ylabel('$\tilde{k}_\mathrm{xx}$ [mD]', latx{:}, 'fontSize', 13)
               xticks([5 10 15 20])
               yticks([1e-3 1e-2 0.1 1 50])
               ylim([7e-4 200])
               grid on
               %legend(p1, {'$\hat{k}_\mathrm{xx}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southeast')
               sgtitle('Med($k$) vs clay residual friction angle', latx{:})
               set(gca,'YScale','log')
               
               subplot(1,3,2)
               hold off
               %xlabel('$\phi^`_\mathrm{r, c}$', 'Interpreter', 'latex', 'fontSize', 13)
               ylabel('$\tilde{k}_\mathrm{yy}$ [mD]', latx{:}, 'fontSize', 13)
               xticks([5 10 15 20])
               yticks([1e-3 1e-2 0.1 1 50])
               ylim([7e-4 200])
               grid on
               set(gca,'YScale','log')
               
               subplot(1,3,3)
               hold off
               %xlabel('$\phi^`_\mathrm{r, c}$', 'Interpreter', 'latex', 'fontSize', 13)
               ylabel('$\tilde{k}_\mathrm{zz}$ [mD]', latx{:}, 'fontSize', 13)
               xticks([5 10 15 20])
               yticks([1e-3 1e-2 0.1 1 50])
               ylim([7e-4 200])
               grid on
               set(gca,'YScale','log')
               
               set(f9, 'position', [500, 500, 900, 300]);
           end
   
       case 'SSFc'
               ssfc{n} = data.SSFc;
           if n == 1
               latx = {'Interpreter','latex'};
               colrs = winter(numel(fname));
               
               f9 = figure(99);
               subplot(1,3,1)
               text(1.1, 1.8e-4, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
               hold on
               text(1.1, 1.8e-3, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
               text(4, 30,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
               loglog([1 12], [data.claySandPerms(1) data.claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               loglog([1 12], [data.claySandPerms(2) data.claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               loglog([1 12], [data.claySandPerms(3) data.claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               
               subplot(1,3,2)
               %text(5, 1.5e-4, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
               hold on
               %text(5, 1.5e-3, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
               %text(5, 30,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
               loglog([1 12], [data.claySandPerms(1) data.claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               loglog([1 12], [data.claySandPerms(2) data.claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               loglog([1 12], [data.claySandPerms(3) data.claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               
               subplot(1,3,3)
               %text(5, 1.5e-4, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
               hold on
               %text(5, 1.5e-3, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
               %text(5, 30,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
               loglog([1 12], [data.claySandPerms(1) data.claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               loglog([1 12], [data.claySandPerms(2) data.claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
               loglog([1 12], [data.claySandPerms(3) data.claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)              
           end
           
           subplot(1,3,1)
           p1 = loglog(ssfc{n}, permMed{n}(:,1), '.-', 'color', colrs(n, :));
           subplot(1,3,2)
           p2 = loglog(ssfc{n}, permMed{n}(:,2), '.-', 'color', colrs(n, :));
           subplot(1,3,3)
           p3 = loglog(ssfc{n}, permMed{n}(:,3), '.-', 'color', colrs(n, :));
           
           if n == numel(fname)
               subplot(1,3,1)
               hold off
               xlabel('$\mathrm{SSF}_\mathrm{c}$', 'Interpreter', 'latex', 'fontSize', 13)
               ylabel('$\tilde{k}_\mathrm{xx}$ [mD]', latx{:}, 'fontSize', 13)
               xlim([1 12])
               xticks([1 2 3 4 5 7 10])
               yticks([1e-4 1e-3 1e-2 0.1 1 10 50])
               yticklabels({'10^{-4}', '10^{-3}', '0.01', '0.1', '1', '10', '50'})
               ylim([7e-5 50])
               grid on
               %legend(p1, {'$\hat{k}_\mathrm{xx}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southeast')
               sgtitle('Med($k$) vs clay residual friction angle', latx{:})
               set(gca,'YScale','log')
               set(gca,'XScale','log')
               
               subplot(1,3,2)
               hold off
               %xlabel('$\phi^`_\mathrm{r, c}$', 'Interpreter', 'latex', 'fontSize', 13)
               ylabel('$\tilde{k}_\mathrm{yy}$ [mD]', latx{:}, 'fontSize', 13)
               xlim([1 12])
               xticks([1 2 3 4 5 7 10])
               yticks([1e-4 1e-3 1e-2 0.1 1 10 50])
               yticklabels({'10^{-4}', '10^{-3}', '0.01', '0.1', '1', '10', '50'})
               ylim([7e-5 50])
               grid on
               set(gca,'YScale','log')
               set(gca,'XScale','log')
               
               subplot(1,3,3)
               hold off
               %xlabel('$\phi^`_\mathrm{r, c}$', 'Interpreter', 'latex', 'fontSize', 13)
               ylabel('$\tilde{k}_\mathrm{zz}$ [mD]', latx{:}, 'fontSize', 13)
               xticks([1 2 3 4 5 7 10])
               yticks([1e-4 1e-3 1e-2 0.1 1 10 50])
               yticklabels({'10^{-4}', '10^{-3}', '0.01', '0.1', '1', '10', '50'})
               ylim([7e-5 50])
               grid on
               xlim([1 12])
               set(gca,'YScale','log')
               set(gca,'XScale','log')
               
               set(f9, 'position', [500, 500, 900, 300]);
           end           
   end
   clear data
end


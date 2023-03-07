function singleParameterPlots(toPlot, nSimSmear, claySandPerms, Perm, Poro, ...
                              clayResidual, sandResidual, SSFc, len, ...
                              faultDip, faultThickness)
%
% single parameter variation analysis
%
latx = {'Interpreter', 'latex'};
switch toPlot        
    case 'permDistr'
        % Permeability distribution
        %idPlot = randi(numel(Perm));
        nbins = 50;
        logMinP = log10(claySandPerms(2)/(milli*darcy));
        logMedP = log10(claySandPerms(3)/(milli*darcy));
        logMaxP = log10(claySandPerms(1)/(milli*darcy));
        edges = logspace(logMinP, logMaxP, nbins);
        
        for n=1:numel(Perm)
            %binProb = histcounts(Perm(n).all(:,1), edges,'Normalization', 'probability');
            %centers = edges(1:end-1)+(diff(edges)/2);
            %[~, permVals] = findpeaks(binProb, centers, 'MinPeakDistance', 1);
            %disp(['Permeability modes: ' num2str(permVals)]);
        
            fh1 = figure(81);
            subplot(3,4,n)
            plot([claySandPerms(1) claySandPerms(1)]./(milli*darcy), [0 1], '-k')
            hold on
            plot([claySandPerms(2) claySandPerms(2)]./(milli*darcy), [0 1], '-.k')
            plot([claySandPerms(3) claySandPerms(3)]./(milli*darcy), [0 1], '--k')
            histogram(Perm(n).all(:,1), edges, 'Normalization', 'probability','FaceColor', [0.3 0.3 0.3])
            ax = gca; ax.XAxis.FontSize = 8;
            if n == 1
                xlabel('$\hat{k}_{xx}$ [mD]', latx{:}, 'fontSize', 13)
                ylabel('P [-]', latx{:}, 'fontSize', 13)
                sgtitle(['Cross-fault permeability distributions'], latx{:}, 'fontSize', 14)
                text(1.1*10^logMinP, 0.9, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
               text(0.1*10^logMedP, 0.2, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
               text(0.3*10^logMaxP, 0.9,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12) 
            end
%             if n == 2
%                text(1.1*10^logMinP, 0.7, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
%                text(1.1*10^logMedP, 0.7, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
%                text(0.3*10^logMaxP, 0.7,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12) 
%             end
            title(['N = ' num2str(nSimSmear(n))])
            hold off
            xlim([0.7*10^logMinP 2*10^logMaxP])
            ylim([0 1])
            grid on
            set(gca,'XScale','log')
            xticks([1e-5 1e-4 1e-3 1e-2 0.1 1 10])
            xticklabels({'10^{-5}' '10^{-4}' '10^{-3}' '0.01' '0.1' '1' '10'})
            set(fh1, 'position', [500, 200, 950, 650]);
            
            
            fh2 = figure(82);
            subplot(3,4,n)
            plot([claySandPerms(1) claySandPerms(1)]./(milli*darcy), [0 1], '-k')
            hold on
            plot([claySandPerms(2) claySandPerms(2)]./(milli*darcy), [0 1], '-.k')
            plot([claySandPerms(3) claySandPerms(3)]./(milli*darcy), [0 1], '--k')
            histogram(Perm(n).all(:,3), edges, 'Normalization', 'probability','FaceColor', [0.3 0.3 1])
            ax = gca; ax.XAxis.FontSize = 8;
            if n == 1
                text(1.1*10^logMinP, 0.9, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
                text(1.1*10^logMedP, 0.9, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
                text(0.3*10^logMaxP, 0.9,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
                xlabel('$\hat{k}_{zz}$ [mD]', latx{:}, 'fontSize', 13)
                ylabel('P [-]', latx{:}, 'fontSize', 13)
                sgtitle(['Updip permeability distributions'], latx{:}, 'fontSize', 14)
            end
            title(['N = ' num2str(nSimSmear(n))])
            hold off
            xlim([0.7*10^logMinP 2*10^logMaxP])
            ylim([0 1])
            grid on
            set(gca,'XScale','log')
            xticks([1e-4 1e-3 1e-2 0.1 1 10])
            xticklabels({'10^{-4}' '10^{-3}' '0.01' '0.1' '1' '10'})
            set(fh2, 'position', [500, 200, 950, 650]);
        end
%         subplot(3,1,2)
%         histogram(Perm(idPlot).all(:,2), edges, 'Normalization', 'probability','FaceColor', [1 0.3 0.3])
%         xlabel('$\hat{k}_{yy}$ [mD]', latx{:}, 'fontSize', 12)
%         ylabel('P [-]', latx{:}, 'fontSize', 13)
%         xlim([10^logMinP 10^logMaxP])
%         ylim([0 1])
%         grid on
%         set(gca,'XScale','log')
%         xticks([1e-3 1e-2 0.1 1 10 100])
    
    case 'numberOfSim'
        % Required N of simulations to assess individual parameter variations
        permMedian = cell2mat(Perm.estimate);
        PermPctErr = cell2mat(Perm.pctErr);
        %Poro = cell2mat(Poro);
        %PoroPctErr = cell2mat(PoroPctErr);
        fh = figure(9);
        loglog([nSimSmear(1) nSimSmear(end)], [claySandPerms(1) claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
        hold on
        loglog([nSimSmear(1) nSimSmear(end)], [claySandPerms(2) claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
        loglog([nSimSmear(1) nSimSmear(end)], [claySandPerms(3) claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
        p1 = errorbar(nSimSmear, permMedian(:,1), PermPctErr(:,1), PermPctErr(:,4), '-ok', 'MarkerSize', 4);
        p2 = errorbar(nSimSmear, permMedian(:,2), PermPctErr(:,2), PermPctErr(:,5), '-or', 'MarkerSize', 4);
        p3 = errorbar(nSimSmear, permMedian(:,3), PermPctErr(:,3), PermPctErr(:,6), '-ob', 'MarkerSize', 4);
        plot(nSimSmear(8), permMedian(8,1), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6)
        plot(nSimSmear(8), permMedian(8,2), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 6)
        plot(nSimSmear(8), permMedian(8,3), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 6)
        hold off
        xlabel('N$_\mathrm{sim}$', 'Interpreter', 'latex', 'fontSize', 14)
        ylabel('$k$ [mD]', latx{:}, 'fontSize', 14)
        xlim([nSimSmear(1) nSimSmear(end)])
        ylim([7e-4 3*10^2])
        grid on 
        legend([p1, p2, p3], {'$\hat{k}_\mathrm{xx}$', '$\hat{k}_\mathrm{yy}$', '$\hat{k}_\mathrm{zz}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'northeast')
        text(25, 1.7e-3, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
        text(25, 7e-3, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
        text(100, 150,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
        title('Median $k$ and [10 90]pp vs N simulations', latx{:})
        set(fh, 'position', [500, 500, 380, 300]);
        set(gca,'YScale','log')
        set(gca,'XScale','log')
    
    case 'clayResidual'
        permMedian = reshape([Perm.estimate], 3, numel(Perm))';
        PermPctErr = reshape([Perm.pctErr], 6, numel(Perm))';
        % Clay residual friction angle
        f10 = figure(10);
        text(14, 1.5e-4, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
        hold on
        text(14, 1.5e-3, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
        text(12, 30,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
        semilogy([clayResidual(1) clayResidual(end)], [claySandPerms(1) claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
        semilogy([clayResidual(1) clayResidual(end)], [claySandPerms(2) claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
        semilogy([clayResidual(1) clayResidual(end)], [claySandPerms(3) claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 1.5)
        p1 = errorbar(clayResidual, permMedian(:,1), PermPctErr(:,1), PermPctErr(:,4), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 4);
        p2 = errorbar(clayResidual, permMedian(:,2), PermPctErr(:,2), PermPctErr(:,5), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 4);
        p3 = errorbar(clayResidual, permMedian(:,3), PermPctErr(:,3), PermPctErr(:,6), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 4);
        hold off
        xlabel('$\phi^`_\mathrm{r, c}$', 'Interpreter', 'latex', 'fontSize', 14)
        ylabel('$k$ [mD]', latx{:}, 'fontSize', 14)
        xticks([5 10 15 20])
        yticks([1e-4 1e-3 1e-2 0.1 1 50])
        ylim([7e-5 50])
        grid on 
        legend([p1, p2, p3], {'$\hat{k}_\mathrm{xx}$', '$\hat{k}_\mathrm{yy}$', '$\hat{k}_\mathrm{zz}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southeast')
        title('Clay residual friction angle', latx{:})
        set(f10, 'position', [500, 500, 400, 350]);
        set(gca,'YScale','log')
        
    case 'SSFc'
        permMedian = reshape([Perm.estimate], 3, numel(Perm))';
        PermPctErr = reshape([Perm.pctErr], 6, numel(Perm))';
        % SSFc
        f10 = figure(10);
        hold on
        text(1.5, 1.5e-4, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
        text(1.5, 1.5e-3, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
        text(7, 30,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
        loglog([SSFc(1) SSFc(end)], [claySandPerms(1) claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        loglog([SSFc(1) SSFc(end)], [claySandPerms(2) claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        loglog([SSFc(1) SSFc(end)], [claySandPerms(3) claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        p1 = errorbar(SSFc, permMedian(:,1), PermPctErr(:,1), PermPctErr(:,4), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 4);
        p2 = errorbar(SSFc, permMedian(:,2), PermPctErr(:,2), PermPctErr(:,5), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 4);
        p3 = errorbar(SSFc, permMedian(:,3), PermPctErr(:,3), PermPctErr(:,6), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 4);
        hold off
        xlabel('SSF$_\mathrm{c}$', 'Interpreter', 'latex', 'fontSize', 13)
        ylabel('$k$ [mD]', latx{:}, 'fontSize', 13)
        xlim([1 max(SSFc)])
        set(gca,'Xtick',[1 2 3 4 5 7 10]);
        yticks([1e-4 1e-3 1e-2 0.1 1 10 50])
        yticklabels({'10^{-4}', '10^{-3}', '0.01', '0.1', '1', '10', '50'})
        ylim([7e-5 50])
        grid on 
        %legend([p1, p2, p3], {'$\hat{k}_\mathrm{xx}$', '$\hat{k}_\mathrm{yy}$', '$\hat{k}_\mathrm{zz}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southwest')
        title('Critical shale smear factor', latx{:}, 'fontsize', 13)
        set(f10, 'position', [500, 500, 500, 400]);
        set(gca,'YScale','log')
        set(gca,'XScale','log')

    case 'sandResidual'
        permMedian = reshape([Perm.estimate], 3, numel(Perm))';
        PermPctErr = reshape([Perm.pctErr], 6, numel(Perm))';
        % Sand residual friction angle
        f10 = figure(10);
        hold on
        text(25.5, 3e-3, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
        text(25.5, 3e-2, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
        text(30, 150,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
        semilogy([sandResidual(1) sandResidual(end)], [claySandPerms(1) claySandPerms(1)]./(milli*darcy), '-k')
        semilogy([sandResidual(1) sandResidual(end)], [claySandPerms(2) claySandPerms(2)]./(milli*darcy), '-.k')
        semilogy([sandResidual(1) sandResidual(end)], [claySandPerms(3) claySandPerms(3)]./(milli*darcy), '--k')
        p1 = errorbar(sandResidual, permMedian(:,1), PermPctErr(:,1), PermPctErr(:,4), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 4);
        p2 = errorbar(sandResidual, permMedian(:,2), PermPctErr(:,2), PermPctErr(:,5), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 4);
        p3 = errorbar(sandResidual, permMedian(:,3), PermPctErr(:,3), PermPctErr(:,6), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 4);
        hold off
        xlabel('$\phi^`_\mathrm{r, s}$', 'Interpreter', 'latex', 'fontSize', 14)
        ylabel('$k$ [mD]', latx{:}, 'fontSize', 14)
        ylim([1e-3 3*10^2])
        grid on 
        %legend([p1, p2, p3], {'$\hat{k}_\mathrm{xx}$', '$\hat{k}_\mathrm{yy}$', '$\hat{k}_\mathrm{zz}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southwest')
        title('Sand residual friction angle', latx{:})
        set(f10, 'position', [500, 500, 400, 350]);
        set(gca,'YScale','log')

    case 'maxSmearLen'
        permMedian = reshape([Perm.estimate], 3, numel(Perm))';
        PermPctErr = reshape([Perm.pctErr], 6, numel(Perm))';
        % max smear length
        f10 = figure(10);
        hold on
        text(10, 1.7e-4, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
        text(10, 6e-4, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
        text(10, 30,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
        loglog([len{1} 100], [claySandPerms(1) claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        loglog([len{1} 100], [claySandPerms(2) claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        loglog([len{1} 100], [claySandPerms(3) claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        p1 = errorbar([[len{1:end-1}] 86.62], permMedian(:,1), PermPctErr(:,1), PermPctErr(:,4), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 4);
        p2 = errorbar([[len{1:end-1}] 86.62], permMedian(:,2), PermPctErr(:,2), PermPctErr(:,5), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 4);
        p3 = errorbar([[len{1:end-1}] 86.62], permMedian(:,3), PermPctErr(:,3), PermPctErr(:,6), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 4);
        hold off
        xlabel('L [m]', 'Interpreter', 'latex', 'fontSize', 13)
        ylabel('$k$ [mD]', latx{:}, 'fontSize', 13)
        xlim([0 100])
        set(gca,'Xtick',[1 2 3 4 5 10 20 30 40 50 100]);
        yticks([1e-4 1e-3 1e-2 0.1 1 10 50])
        yticklabels({'10^{-4}', '10^{-3}', '0.01', '0.1', '1', '10', '50'})
        ylim([7e-5 50])
        grid on 
        %legend([p1, p2, p3], {'$\hat{k}_\mathrm{xx}$', '$\hat{k}_\mathrm{yy}$', '$\hat{k}_\mathrm{zz}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southwest')
        title('Max smear segment length', latx{:}, 'fontSize', 13)
        set(f10, 'position', [500, 500, 400, 350]);
        set(gca,'YScale','log')
        set(gca,'XScale','log')

    case 'faultDip'
        % Fault dip
        f10 = figure(10);
        hold on
        text(2.5, 3e-3, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
        text(2.5, 3e-2, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
        text(10, 150,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
        semilogy([0 90], [claySandPerms(1) claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        semilogy([0 90], [claySandPerms(2) claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        semilogy([0 90], [claySandPerms(3) claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        p1 = errorbar(faultDip, Perm(:,1), PermPctErr(:,1), PermPctErr(:,4), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 4);
        p2 = errorbar(faultDip, Perm(:,2), PermPctErr(:,2), PermPctErr(:,5), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 4);
        p3 = errorbar(faultDip, Perm(:,3), PermPctErr(:,3), PermPctErr(:,6), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 4);
        hold off
        xlabel('$\beta$ $[^\circ]$', 'Interpreter', 'latex', 'fontSize', 14)
        ylabel('$k$ [mD]', latx{:}, 'fontSize', 14)
        xlim([0 90])
        set(gca,'Xtick',[10 20 30 40 50 60 70 80 90]);
        ylim([1e-3 3*10^2])
        grid on 
        %legend([p1, p2, p3], {'$\hat{k}_\mathrm{xx}$', '$\hat{k}_\mathrm{yy}$', '$\hat{k}_\mathrm{zz}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southwest')
        title('Fault dip', latx{:})
        set(f10, 'position', [500, 500, 400, 350]);
        set(gca,'YScale','log')

    case 'faultThick'
        permMedian = reshape([Perm.estimate], 3, numel(Perm))';
        PermPctErr = reshape([Perm.pctErr], 6, numel(Perm))';
        % Fault thickness
        f10 = figure(10);
        hold on
        text(4, 1.7e-4, '$k_\mathrm{c,\perp}$', latx{:}, 'fontSize', 12)
        text(4, 1.7e-3, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
        text(0.3, 30,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
        loglog([0.1 10], [claySandPerms(1) claySandPerms(1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        loglog([0.1 10], [claySandPerms(2) claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        loglog([0.1 10], [claySandPerms(3) claySandPerms(3)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        p1 = errorbar(faultThickness, permMedian(:,1), PermPctErr(:,1), PermPctErr(:,4), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 4);
        p2 = errorbar(faultThickness, permMedian(:,2), PermPctErr(:,2), PermPctErr(:,5), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 4);
        p3 = errorbar(faultThickness, permMedian(:,3), PermPctErr(:,3), PermPctErr(:,6), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 4);
        hold off
        xlabel('f$_\mathrm{T}$ [m]', 'Interpreter', 'latex', 'fontSize', 13)
        ylabel('$k$ [mD]', latx{:}, 'fontSize', 13)
        xlim([0.1 10])
        set(gca,'Xtick',[0.1 0.5 1 1.5 2 3 4 5 7 10]);
        yticks([1e-4 1e-3 1e-2 0.1 1 10 50])
        yticklabels({'10^{-4}', '10^{-3}', '0.01', '0.1', '1', '10', '50'})
        ylim([7e-5 50])
        grid on 
        %legend([p1, p2, p3], {'$\hat{k}_\mathrm{xx}$', '$\hat{k}_\mathrm{yy}$', '$\hat{k}_\mathrm{zz}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southwest')
        title('Fault thickness', latx{:}, 'fontSize', 13)
        set(f10, 'position', [500, 500, 400, 350]);
        set(gca,'YScale','log')
        set(gca,'XScale','log')

    case 'clayPerms'
        permMedian = reshape([Perm.estimate], 3, numel(Perm))';
        %PermPctErr = reshape([Perm.pctErr], 6, numel(Perm))';
        cp = claySandPerms(:,2:3)/(milli*darcy);
        % ClayPerms
        f10 = figure(10);
        hold on
        text(2, 2e-5, 'min$(k_\mathrm{c,\perp})$', latx{:}, 'fontSize', 12)
        %text(2, 1e-4, '$k_\mathrm{c,\backslash\backslash}$', latx{:}, 'fontSize', 12)
        text(2, 30,  '$k_\mathrm{s}$', latx{:}, 'fontSize', 12)
        semilogy([1 4], [claySandPerms(1,1) claySandPerms(1,1)]./(milli*darcy), '-', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        %semilogy([1 5], [claySandPerms(2) claySandPerms(2)]./(milli*darcy), '-.', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        semilogy([1 4], [claySandPerms(end,2) claySandPerms(end,2)]./(milli*darcy), '--', 'color', [0.5 0.5 0.5], 'linewidth', 2)
        semilogy(1:4, permMedian(:,1)./cp(:,1), '-ok', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 4);
        semilogy(1:4, permMedian(:,2)./cp(:,2), '-or', 'MarkerFaceColor', [255, 125, 125]/255, 'MarkerSize', 4);
        semilogy(1:4, permMedian(:,3)./cp(:,2), '-ob', 'MarkerFaceColor', [125, 125, 255]/255, 'MarkerSize', 4);
        hold off
        xlabel('$k_\mathrm{c}$ [mD]', 'Interpreter', 'latex', 'fontSize', 14)
        ylabel('$\kappa$', latx{:}, 'fontSize', 14)
        xlim([1 4])
        set(gca,'Xtick',1:4);
        xticklabels({'$\frac{k_\mathrm{\perp}=0.01}{k_\mathrm{\backslash\backslash}=0.1}$', ...
                     '$\frac{10^{-3}}{10^{-2}}$','$\frac{10^{-4}}{10^{-3}}$',...
                     '$\frac{10^{-5}}{10^{-4}}$'})
        set(gca,'TickLabelInterpreter','latex')
        ylim([1 1e5])
        %yticks([1e-5 1e-4 1e-3 1e-2 0.1 1 50])
        grid on 
        %legend([p1, p2, p3], {'$\hat{k}_\mathrm{xx}$', '$\hat{k}_\mathrm{yy}$', '$\hat{k}_\mathrm{zz}$'}, latx{:}, 'fontSize', 12, 'box', 'on', 'location', 'southwest')
        title('Normalized output permeability', latx{:})
        set(f10, 'position', [500, 500, 400, 350]);
        set(gca,'YScale','log')
        set(gca,'XDir','reverse')
end

end
function [f, rock, cP, fault, refPc] = getFluid_bo(mesh, G, rock, opt, ucids, fault, figs)
%
%
%

% Assign saturation and cp regions
if strcmp(mesh.type, 'sc2') && mesh.ampLyr == 21
    if opt.sc == 1
        id_ampb = 2:22;
    elseif opt.sc == 2
        id_ampb = 2:2:22;
    end
    id_ampb  = [id_ampb 58];
    id_seal  = [ucids{id_ampb}];
    %id_res   = [ucids{1, 3:4}];
    id_fault = ucids{61};
    id_faultRock = [ucids{27:30}];
    fault.satRegNum = 3;
    res_satRegNum = 1;
    
% elseif strcmp(mesh.type, 'coarse') || strcmp(mesh.type, 'ref300') || ...
%    strcmp(mesh.type, 'ref200') || strcmp(mesh.type, 'ref100') || ...
%     strcmp(mesh.type, 'ref50') || strcmp(mesh.type, 'ref40')
%     id_seal  = [ucids{1} ucids{3} ucids{5} ucids{9} ...
%                 ucids{10} ucids{11}];
%     %id_res   = [ucids{2} ucids{4} ucids{6} ucids{8} ...
%     %            ucids{12} ucids{13}];
%     id_fault = [];
%     id_faultRock = [];
%     fault.satRegNum = 2;
%     id_ampb = 5;
% elseif strcmp(mesh.type, 'ref25') && mesh.reduce == 1 || ...
%        strcmp(mesh.type, 'ref20') && mesh.reduce == 1
%     id_seal  = [ucids{2} ucids{6} ucids{7} ucids{8}];
%     id_fault = [];
%     id_faultRock = [];
%     fault.satRegNum = 2;
%     id_ampb = 2;    
% elseif strcmp(mesh.type, 'ref15') && mesh.reduce == 2
%     id_seal  = [ucids{2} ucids{4} ucids{5}];
%     id_fault = [];
%     id_faultRock = [];
%     fault.satRegNum = 2;
%     id_ampb = 2;   
% else % upper
%         id_seal  = ucids{2};
%         %id_res   = [ucids{1, 3:4}];
%         id_fault = ucids{end};
%         id_faultRock = [ucids{6:8}];
%         fault.satRegNum = 3;
%        id_ampb = 2;
end
seal_satRegNum = 2;
% 3 regions for saturation, to account for different capillary properties
% in the fault
rock.regions.saturation = ones(G.cells.num,1);                              % reservoir units
rock.regions.saturation(id_seal) = seal_satRegNum;
if strcmp(opt.faultPc, 'section_predict')
    % MESH DEPENDENT
    rock.regions.saturation(ucids{26}) = 3;
    rock.regions.saturation(ucids{27}) = 4;
    rock.regions.saturation(ucids{28}) = 5;
    rock.regions.saturation(ucids{29}) = 6;
    rock.regions.saturation(ucids{30}) = 7;
    rock.regions.saturation(ucids{31}) = 8;
else
    rock.regions.saturation(id_fault) = fault.satRegNum;
end

% 2 regions for compressibility
rock.regions.rocknum = ones(G.cells.num,1); 
rock.regions.rocknum([id_seal, id_faultRock], 1) = 2;                       % only required for keyword ROCK

% Read data and define fluid object
switch opt.fluid
    case 'sc2_krHyst_revised'           % SC2 BaseCase, DEFSIM
        fn   = 'ls-proj/gom/adbo_2.5D/eclipse_data_files/sc2_krHyst_fPcFault_revised.DATA';
    case 'krHyst_pcFaultScaled_pvtBase' % OLD
        %fn   = 'ls-proj/gom/josimar_coupled_3D/data/josimarGoM_hyst.DATA';  % used for Josimar's paper plots!
        fn   = 'ls-proj/gom/adbo_2.5D/eclipse_data_files/krHyst_fPcFault_co2brine.DATA';
    case 'sc1_krHyst_revised'           % SC1, DEFSIM
        fn = 'ls-proj/gom/adbo_2.5D/eclipse_data_files/sc1_krHyst_revised.DATA';
    case 'krHyst_FaultPcKrUps'          % PREDICT-based upscaling
        if strcmp(opt.permCase, 'high')
            fn   = 'ls-proj/gom/adbo_2.5D/eclipse_data_files/krHyst_AmpBSectionRegions_highPerm.DATA';
        else
            fn   = 'ls-proj/gom/adbo_2.5D/eclipse_data_files/krHyst_AmpBSectionRegions_lowPerm.DATA'; 
        end
    case 'noHyst_FaultPcKrUps'
        if strcmp(opt.permCase, 'high')
            fn   = 'ls-proj/gom/adbo_2.5D/eclipse_data_files/noHyst_AmpBSectionRegions_highPerm.DATA';
        else
            fn   = 'ls-proj/gom/adbo_2.5D/eclipse_data_files/noHyst_AmpBSectionRegions_lowPerm.DATA'; 
        end
end
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.ROCKNUM = rock.regions.rocknum;
f = initDeckADIFluid(deck);

% Get pore compressibility
cP = deck.PROPS.ROCK(:, 2);

% Scaled Pc in fault
if strcmp(opt.faultPc, 'scaled')
   %refPc.val  = f.pcOG{seal_satRegNum};
   refPc.val  = f.pcOG{fault.satRegNum};
   refPc.poro = 0.07;                            % Table 3.3, sample 5 in Treviño & Meckel, 2017 (Pc curve is the green one in Fig. 3.4a)
   refPc.perm = 0.00208*milli*darcy;             % "
else
   refPc = [];
end

% Hysteresis
if isfield(f, 'krHyst') && f.krHyst == 1
   numReg = max(rock.regions.saturation);
   f.krHyst = res_satRegNum + numReg;            % porous rock is first imb region (4), fault is third (6)
   rock.regions.imbibition = rock.regions.saturation + numReg;
   f = addScanKr(f, rock.regions.imbibition, 0.05);
end

% Figures
if figs == 1
    plotFluid(f)
end


%% Helper functions
    function plotFluid(f)
        latx  = {'Interpreter','latex'};
        clrs  = [65,2,0; 164,30,53; 125,52,162; 49,148,206; 220,220,220]./255;
        Sbr   = linspace(f.krPts.og(1, 2), f.krPts.og(1, 3), 50);
        Sbs   = linspace(f.krPts.og(2, 2), f.krPts.og(2, 3), 50);
        % background color for print: [253 253 248]/255
        
        % kr, Pc Reservoir
        h = figure(30);
        plot(Sbr, f.krOG{1}(Sbr), '-', 'color', clrs(4,:), 'linewidth', 1.5);
        hold on; plot(Sbr, f.krG{1}(1-Sbr), '-', 'color', clrs(3,:), 'linewidth', 1.5); 
        plot(Sbr, f.krG{4}(1-Sbr), '--', 'color', clrs(3,:), 'linewidth', 1.5);

        
        hold off, xlim([0 1]), grid on;
        ax = gca; ax.YAxis.FontSize = 12; ax.XAxis.FontSize = 12;
        xlabel('S$_\mathrm{aq}$ [-]', latx{:}, 'fontsize', 15)
        ylabel('k$_\mathrm{r}$ [-]','fontsize', 15, latx{:})
        ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
        set(ax,'Xtick',[0 0.2 0.4 0.6 0.8 1])
        set(ax,'Ytick',[0 0.25 0.5 0.75 1])
        legend('k$_\mathrm{r,aq}$', 'k$_\mathrm{r,g}$', 'k$_\mathrm{r,g}^\mathrm{imb}$', ...
               latx{:}, 'location', 'north', 'box', 'off', ...
               'fontsize', 14, 'textcolor', [100 100 100]/255)
        set(h, 'Position', [600, 600, 380, 320])
        %set(ax,'color','none')
        
        
        % kr, Pc Seal
        h = figure(31);
        yyaxis left
        plot(Sbs, f.krOG{2}(Sbs), '-', 'color', clrs(4,:), 'linewidth', 1.5);
        hold on; plot(Sbs, f.krG{2}(1-Sbs), '-', 'color', clrs(3,:), 'linewidth', 1.5); 
        hold off,  xlim([0 1]), grid on; ax = gca;
        xlabel('S$_\mathrm{aq}$ [-]','fontsize', 15, latx{:})
        ylabel('k$_\mathrm{r}$ [-]','fontsize', 15, latx{:})
        set(ax,'Xtick',[0 0.2 0.4 0.6 0.8 1]); set(ax,'Ytick',[0 0.25 0.5 0.75 1])
        yyaxis right
        plot(Sbs, f.pcOG{2}(1-Sbs)/barsa, 'color', clrs(3,:), 'LineStyle', '-.', 'linewidth', 1.5); 
        ylabel('P$_\mathrm{c}$ [bar]','fontsize', 15, latx{:})
        ax.YAxis(1).Color = [100 100 100]/255; ax.YAxis(1).FontSize = 12;
        ax.YAxis(2).Color = [100 100 100]/255; ax.YAxis(2).FontSize = 12;
        ax.XAxis.Color = [100 100 100]/255; ax.XAxis.FontSize = 12;
        legend('k$_\mathrm{r,aq}$', 'k$_\mathrm{r,g}$', 'P$_\mathrm{c}$', ...
               latx{:}, 'location', 'northwest', 'box', 'off', ...
               'fontsize', 14, 'textcolor', [100 100 100]/255)
        set(h, 'Position', [600, 600, 380, 320])
        
        % Fault
        if numel(f.krOG) == 3
            Sbf   = linspace(f.krPts.og(3, 2), f.krPts.og(3, 3), 50);
            h = figure(32);
            plot(Sbf, f.krOG{3}(Sbf), '-', 'color', clrs(4,:), 'linewidth', 1.5);
            hold on; plot(Sbf, f.krG{3}(1-Sbf), '-', 'color', clrs(3,:), 'linewidth', 1.5); 
            hold off, xlim([0 1]), grid on;
            ax = gca; ax.YAxis.FontSize = 12; ax.XAxis.FontSize = 12;
            xlabel('S$_\mathrm{aq}$ [-]', latx{:}, 'fontsize', 15)
            ylabel('k$_\mathrm{r}$ [-]','fontsize', 15, latx{:})
            ax.YAxis.Color = [100 100 100]/255; ax.XAxis.Color = [100 100 100]/255;
            set(ax,'Xtick',[0 0.2 0.4 0.6 0.8 1])
            set(ax,'Ytick',[0 0.25 0.5 0.75 1])
            legend('k$_\mathrm{r,aq}$', 'k$_\mathrm{r,g}$', ...
                   latx{:}, 'location', 'north', 'box', 'off', ...
                   'fontsize', 14, 'textcolor', [100 100 100]/255)
            set(h, 'Position', [600, 600, 380, 320])
        end
        
        % Density, viscosity
        p = linspace(1, 800, 50)*barsa;
        rhoG = f.bG(p).*f.rhoGS;
        muG  = f.muG(p);
        rs   = f.rsSat(p);
        rhoB_sat = f.bO(p, rs, true(numel(p), 1)).*(rs.*f.rhoGS + f.rhoOS);
        rhoB = f.bO(p',zeros(50,1),false(50,1))*f.rhoOS;
        muB  = f.muO(p, rs, true(numel(p), 1));
        
        col_co2  = 'r'; %[3, 130, 30]./255;
        col_brine  = 'b'; %[0, 71, 148]./255;
        
        h = figure(35);
        yyaxis left
        %plot(p/mega, rhoB_sat, '-', 'color', col_brine, 'linewidth', 1.5);
        hold on;
        plot(p/mega, rhoB, '-', 'color', col_brine, 'linewidth', 1);
        plot(p/mega, rhoG, '-', 'color', col_co2, 'linewidth', 1); 
        hold off,  xlim([0 max(p)/barsa]), grid on; ax = gca;
        xlabel('$p$ [MPa]','fontsize', 14, latx{:})
        ylabel('$\rho$ [kg/m$^3$]','fontsize', 14, latx{:})
        set(ax,'Xtick',[0 10 20 30 40 50]); xlim([0 50])
        set(ax,'Ytick',[0 250 500 750 1000 1100]); ylim([0 1100]);
        yyaxis right
        plot(p/mega, muB, 'color', col_brine, 'LineStyle', '--', 'linewidth', 1);
        hold on; plot(p/mega, muG, '--', 'color', col_co2, 'linewidth', 1); 
        ylabel('$\mu$ [Pa$\cdot$s]','fontsize', 14, latx{:})
        set(ax,'Ytick',[0 0.15 0.3 0.45 0.6 0.66]*10^-3); ylim([0 0.66]*10^-3)
        %ax.YAxis(1).Color = [100 100 100]/255; 
        ax.YAxis(1).Color = 'k';
        ax.YAxis(1).FontSize = 11;
        ax.YAxis(2).Color = 'k'; 
        ax.YAxis(2).FontSize = 11;
        ax.XAxis.Color = 'k'; 
        ax.XAxis.FontSize = 11;
        legend('$\rho_\mathrm{aq}$', '$\rho_\mathrm{g}$', '$\mu_\mathrm{aq}$', '$\mu_\mathrm{g}$', ...
               latx{:}, 'location', 'east', 'box', 'off', ...
               'fontsize', 12) %, 'textcolor', [100 100 100]/255)
        set(h, 'Position', [600, 600, 400, 350])
        
    end

end

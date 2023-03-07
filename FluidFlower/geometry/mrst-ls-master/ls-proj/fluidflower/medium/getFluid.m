function [fluid, rock] = getFluid(inj_type, model_case, G, G_dat, rock, unit, opt, plts, makePlot)
%
%
%
if inj_type == 1
    if model_case == 1
        pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/fluid_props/model1/');
        %fn  = fullfile(pth, 'props_model1.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0.1.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0.01.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.5_Cpc0.5.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.25_Cpc0.5.DATA');   
        fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125_CPc0.5.DATA');  % match model 1
    elseif model_case == 2
        pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/fluid_props/model2/');
        %fn  = fullfile(pth, 'props_model2.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0.1.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.5_Cpc0.5.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.25_Cpc0.5.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.125_Cpc0.33.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.18_Cpc0.33.DATA');
        fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.15_Cpc0.33.DATA');  % match model 2
    elseif model_case == 3
        pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/fluid_props/model3/');
        %fn  = fullfile(pth, 'props_model3.DATA');
        %fn  = fullfile(pth, 'props_model3_Cpc1.5.DATA');
        fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5.DATA'); % match model 3
    end
elseif inj_type == 2
    if model_case == 1
        pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/fluid_props/model1/');
        %fn  = fullfile(pth, 'props_model1.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0.1.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0.01.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.5_Cpc0.5.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.25_Cpc0.5.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125_CPc0.5.DATA');  % match model 1
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125_CPc0.5_CfinfPce7.3mb_CfsupPce4.5mb.DATA'); % AC07
        fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125_CPc0.5_CfinfPce5mb_CfsupPce3.5mb.DATA'); %AC07
    elseif model_case == 2
        pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/fluid_props/model2/');
        %fn  = fullfile(pth, 'props_model2.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0.1.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.5_Cpc0.5.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.25_Cpc0.5.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.125_Cpc0.33.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.18_Cpc0.33.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.15_Cpc0.33.DATA');  % match model 2
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.15_Cpc0.33_CfinfPce7.3mb_CfsupPce4.5mb.DATA'); %AC07
        fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.15_Cpc0.33_CfinfPce5mb_CfsupPce3.5mb.DATA'); %AC07
    elseif model_case == 3
        pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/fluid_props/model3/');
        %fn  = fullfile(pth, 'props_model3.DATA');
        %fn  = fullfile(pth, 'props_model3_Cpc1.5.DATA');
        %fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5.DATA'); % match model 3
        %fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5_CfinfPce7.3mb_CfsupPce4.5mb.DATA'); %AC07
        %fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5_CfinfPce4.5mb_CfsupPce4.5mbSge1e-4.DATA'); %AC07
        %fn  = fullfile(pth, 'props_model3_ESFpc2_CfinfPce5mb_CrestPce3.5mbSge1e-4.DATA'); %AC07
        %fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5_CfinfPce5mb_CfsupPce3.5mb.DATA'); %AC07
        fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5Sge1e-4_CfinfPce5mbSge1e-3_CfsupPce3.5mbSge1e-3.DATA');
        %fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5Sge1e-5_CfinfPce5mb_CfsupPce3.5mb.DATA');
        %fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5Sge1e-5_CfinfPce5mbSge1e-4_CfsupPce3.5mb.DATA');
    end
elseif inj_type == 3
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/medium/fluid_props/Bilbo/');
    if model_case == 1
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125_Cpc0.5.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125.DATA');
        fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125_CPceSg1e-4.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125_CPc2_CPceSg1e-5.DATA');
        %fn  = fullfile(pth, 'props_model1_FPc0_EPc0.125_CPc2_CPceSg1e-4.DATA');
    elseif model_case == 2
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.15_Cpc0.33.DATA');
        fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.15_Cpc1.3.DATA');
        %fn  = fullfile(pth, 'props_model2_Fpc0_Epc0.15_Cpc2.DATA');
    elseif model_case == 3
        %fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc1.5.DATA');
        fn  = fullfile(pth, 'props_model3_ESFpc2_Cpc3.2.DATA');
    end 
end

% Saturation regions
if inj_type == 1
    rock.regions.saturation = ones(G.cells.num, 1);    % Fsup and Finf (and Fmid)
    rock.regions.saturation(ismember(G_dat.p, unit.E)) = 2;
    rock.regions.saturation(ismember(G_dat.p, unit.C)) = 3;
    rock.regions.saturation(ismember(G_dat.p, unit.Cf)) = 3;
    rock.regions.saturation(ismember(G_dat.p, unit.ESF)) = 4;
elseif inj_type == 2
    rock.regions.saturation = ones(G.cells.num, 1);    % Fsup and Finf (and Fmid)
    rock.regions.saturation(ismember(G_dat.p, unit.E)) = 2;
    rock.regions.saturation(ismember(G_dat.p, unit.C)) = 3; % 3
    rock.regions.saturation(ismember(G_dat.p, unit.Cf(1))) = 4;     %Cfsup 4
    rock.regions.saturation(ismember(G_dat.p, unit.Cf(2))) = 5;     %Cfinf 5
    rock.regions.saturation(ismember(G_dat.p, unit.ESF)) = 6;       % 6
elseif inj_type == 3
    rock.regions.saturation = ones(G.cells.num, 1);    % Fsup and Finf (and Fmid)
    rock.regions.saturation(ismember(G_dat.p, unit.Esup)) = 2;
    rock.regions.saturation(ismember(G_dat.p, unit.Einf)) = 2;
    rock.regions.saturation(ismember(G_dat.p, unit.Csup)) = 3;
    rock.regions.saturation(ismember(G_dat.p, unit.Cinf)) = 3;
    rock.regions.saturation(ismember(G_dat.p, unit.CESF)) = 3;
    rock.regions.saturation(ismember(G_dat.p, unit.ESF)) = 4;
end

% Regions for rock compressibility
rock.regions.rocknum = ones(G.cells.num,1); 

% Load deck and initialize fluid
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.ROCKNUM = rock.regions.rocknum;
fluid = initDeckADIFluid(deck);

% Get pore compressibility
cP = deck.PROPS.ROCK(:, 2);
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);

% Hysteresis
if strcmp(opt.hyster, 'on') && isfield(fluid, 'krHyst')
   numReg = max(rock.regions.saturation);
   fluid.krHyst = [1 2 3] + numReg;                % all but ESF
   rock.regions.imbibition = rock.regions.saturation + numReg;
   fluid = addScanKr(fluid, rock.regions.imbibition, 0.02);
end

if makePlot == 1
    plts;
   % Plot density, viscosity
    % f = fluid;
    % p = linspace(1, 5, 50)*barsa;
    % rhoG = f.bG(p).*f.rhoGS;
    % muG  = f.muG(p);
    % rs   = f.rsSat(p);
    % rhoB = f.bO(p, rs, true(numel(p), 1)).*(rs.*f.rhoGS + f.rhoOS);
    % rhoB_noCO2 = f.bO(p', zeros(numel(p), 1), false(numel(p), 1)).*f.rhoOS;
    % muB  = f.muO(p, rs, true(numel(p), 1));
    % plot(p, rhoB, '-b'); hold on; plot(p, rhoB_noCO2, '--b')
    % 
    % p2 = 4*barsa*ones(50,1);
    % rs = linspace(0,f.rsSat(p2(1)),50)';
    % rhoB2 = [f.bO(p2(1:end-1),rs(1:end-1),false(49,1)); ...
    %          f.bO(p2(1),rs(end),true)].*(rs.*f.rhoGS + f.rhoOS);
    % plot(rs, rhoB2, '-k'); xlabel('rs'), ylabel('\rho [kg/m^3]')

    % Plot kr, Pc
    f = figure(11);
    nreg = max(rock.regions.saturation);
    latx = {'Interpreter', 'latex'};
    fontsz = {'fontSize', 16};
    sand_type = {'F', 'E', 'C', 'ESF'};
    for n=1:nreg
        sw = 1 - [0 1e-3 0.01:0.02:fluid.krPts.g(n,3)];
        kw = fluid.krOG{n}(sw);
        kg = fluid.krG{n}(1-sw);
        kgi = fluid.krG{n+nreg}(1-sw);
        h= subplot(1,nreg,n);
        hold on
        plot(sw,kw, '-b', 'linewidth', 1)
        plot(sw,kg, '-r', 'linewidth', 1)
        plot(sw,kgi, '--r', 'linewidth', 1)
        hold off
        grid on
        if n == 1
            xlabel('$S_w$ [-]', latx{:}, fontsz{:})
            ylabel('$k_r$ [-]', latx{:}, fontsz{:})
            legend({'$k_{r,w}$','$k_{r,g}$','$k_{r,g}^{i}$'}, latx{:}, 'fontSize', 12)
        end
        xticks([0 0.2 0.4 0.6 0.8 1])
        yticks([0 0.2 0.4 0.6 0.8 1])
        title(sand_type{n}, fontsz{:})
        h.FontSize = 14;
    end
    f.Position = [0, 0, 1000, 250];

    f = figure(12);
    for n=2:nreg
        sg = [0 1e-3 0.01:0.01:fluid.krPts.g(n,3)];
        pc = fluid.pcOG{n}(sg)/100;
        h = subplot(1,nreg-1,n-1);
        plot(1-sg, pc, '-k', 'LineWidth', 1)
        grid on
        xlabel('$S_w$ [-]', latx{:}, fontsz{:})
        ylabel('$P_c$ [mbar]', latx{:}, fontsz{:})
        xticks([0 0.2 0.4 0.6 0.8 1])
        title(sand_type{n}, fontsz{:})
        h.FontSize = 14;
    end
    f.Position = [0, 0, 750, 250];
end
end
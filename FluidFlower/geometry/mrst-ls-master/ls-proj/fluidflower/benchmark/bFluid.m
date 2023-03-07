function [fluid, rock] = bFluid(inj_case, model_case, G, G_dat, rock, ...
                                unit, removeS, revisedPc, opt, plts, makePlot)
%
%
%
if opt.temp_C == 20
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/benchmark/fluid_props/20C');
elseif opt.temp_C == 25
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/benchmark/fluid_props/25C');
elseif opt.temp_C == 2525
    pth = fullfile(mrstPath('ls-proj'), 'fluidflower/benchmark/fluid_props/25C_25std');
end
if model_case == 1
    if removeS
        %fn  = fullfile(pth, 'bprops_model1_ESFPce1e-3_Sremoved.DATA');
        if opt.pcD == 0.3
            if opt.modelParams == 1
                fn  = fullfile(pth, 'bprops_model1_Sremoved_pcD0.3_1.DATA');
            elseif opt.modelParams == 2
                fn  = fullfile(pth, 'bprops_model1_Sremoved_pcD0.3_2.DATA');
            elseif opt.modelParams == 3
                fn  = fullfile(pth, 'bprops_model1_Sremoved_pcD0.3.DATA');
            end
        else
            fn  = fullfile(pth, 'bprops_model1_Sremoved.DATA');
        end
    else
        fn  = fullfile(pth, 'bprops_model1.DATA');
    end
    
elseif model_case == 2
    if removeS
        %fn  = fullfile(pth, 'bprops_model2_Sremoved.DATA');
        %fn  = fullfile(pth, 'bprops_model22_Sremoved.DATA');    % modif krg,w
        if opt.modelParams == 1
            fn  = fullfile(pth, 'bprops_model2_Sremoved_1.DATA'); 
        elseif opt.modelParams == 2
            fn  = fullfile(pth, 'bprops_model2_Sremoved_2.DATA'); 
        elseif opt.modelParams == 3
            fn  = fullfile(pth, 'bprops_model2_Sremoved.DATA');   
        end
    else
        fn  = fullfile(pth, 'bprops_model2.DATA');
    end
elseif model_case == 3
    if removeS
        if revisedPc == 0
            if opt.modelParams == 1
                fn  = fullfile(pth, 'bprops_model3_Sremoved_1.DATA'); 
            elseif opt.modelParams == 2
                fn  = fullfile(pth, 'bprops_model3_Sremoved_2.DATA'); 
            elseif opt.modelParams == 3
                fn  = fullfile(pth, 'bprops_model3_Sremoved.DATA');    
            end
        elseif revisedPc == 1
            fn  = fullfile(pth, 'bprops_model3_Sremoved_adjustPc_CDE.DATA');
        elseif revisedPc == 2
            fn  = fullfile(pth, 'bprops_model3_Sremoved_adjustPc_CE.DATA');
        end
    end
end

% Saturation regions
rock.regions.saturation = ones(G.cells.num, 1);            % G
rock.regions.saturation(ismember(G_dat.p, unit.F)) = 2;    % F
rock.regions.saturation(ismember(G_dat.p, unit.E)) = 3;    % E
rock.regions.saturation(ismember(G_dat.p, unit.D)) = 4;    % D
rock.regions.saturation(ismember(G_dat.p, unit.C)) = 5;    % C
rock.regions.saturation(ismember(G_dat.p, unit.ESF)) = 6;  % ESF
rock.regions.saturation(ismember(G_dat.p, unit.ESFsup)) = 6;  % ESF
if ~removeS
    rock.regions.saturation(ismember(G_dat.p, unit.S)) = 7;    % S
end

% Regions for rock compressibility
rock.regions.rocknum = ones(G.cells.num,1); 

% Load deck and initialize fluid
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.ROCKNUM = rock.regions.rocknum;
fluid = initDeckADIFluid(deck);


% Get pore compressibility
cP = deck.PROPS.ROCK(:, 2);
%cP = 0;
fluid = assignPvMult(fluid, cP, rock.regions.rocknum);

% Hysteresis
if strcmp(opt.hyster, 'on') && isfield(fluid, 'krHyst')
   disp('****************************************************************')
   disp('********************* HYSTERESIS ACTIVE ************************')
   disp('****************************************************************')
   numReg = max(rock.regions.saturation);
   fluid.krHyst = [1 2 3 4 5] + numReg;                % all but ESF
   rock.regions.imbibition = rock.regions.saturation + numReg;
   fluid = addScanKr(fluid, rock.regions.imbibition, opt.minSat);
else  % remove hyst fields for welltest
   fluid = rmfield(fluid, 'krHyst');
   fluid = rmfield(fluid, 'ehystr');
end

if makePlot == 1
    %plts;
    
    % Plot density, viscosity
    f = fluid;
    p = linspace(1, 5, 50)*barsa;
    rhoG = f.bG(p).*f.rhoGS;
    muG  = f.muG(p);
    rs   = f.rsSat(p);
    rhoB = f.bO(p, rs, true(numel(p), 1)).*(rs.*f.rhoGS + f.rhoOS);
    rhoB_noCO2 = f.bO(p', zeros(numel(p), 1), false(numel(p), 1)).*f.rhoOS;
    muB  = f.muO(p, rs, true(numel(p), 1));
    plot(p, rhoB, '-b'); hold on; plot(p, rhoB_noCO2, '--b')
    
    p2 = 4*barsa*ones(50,1);
    rs = linspace(0,f.rsSat(p2(1)),50)';
    rhoB2 = [f.bO(p2(1:end-1),rs(1:end-1),false(49,1)); ...
             f.bO(p2(1),rs(end),true)].*(rs.*f.rhoGS + f.rhoOS);
    plot(rs, rhoB2, '-k'); xlabel('rs'), ylabel('\rho [kg/m^3]')

    % Plot kr, Pc
    
end
end
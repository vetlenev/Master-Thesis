function rock = bRock(model_case, G, G_dat, unit, m, removeS, opt, plts, plotPoro, plotPerm)
%
%
%

% Diameters are just for model_case 1 (no access to poro-perm data)
%       S,   ESF,  C,    D,    E,   F,   G
%d_mm = [nan, 0.1, 0.9, 1.15, 1.7, opt.dF, 3.25]; % average sand-grain diameters [mm]
if opt.modelParams == 1
    d_mm = [nan, 0.2, 0.66, 1.05, 1.45, 1.77, 2.51];
elseif opt.modelParams == 2
    d_mm = [nan, 0.2, 0.66, 1.05, 1.45, 1.77, 2.51];
elseif opt.modelParams == 3
    d_mm = [nan, 0.1, 0.66, 1.15, 1.45, opt.dF, 3.25];
end
if model_case == 1      % Best match medium fluidflower
    % --------------------------------------------------------------------
    % Porosity
    % --------------------------------------------------------------------
    % Data from Beard & Weyl, AAPG (1973) for moderately to well
    % sorted and from Smits et al. (2010).
    porov = [0.01, .37, .38, .41, .39, .39, .42];    % S, ESF, C, D, E, F, G 
    rock.poro = nan(G.cells.num, 1);
    if ~removeS
        rock.poro(ismember(G_dat.p, unit.S)) = porov(1);
    end
    rock.poro(ismember(G_dat.p, unit.ESF)) = porov(2);
    rock.poro(ismember(G_dat.p, unit.ESFsup)) = porov(2);
    rock.poro(ismember(G_dat.p, unit.C)) = porov(3);
    rock.poro(ismember(G_dat.p, unit.D)) = porov(4);
    rock.poro(ismember(G_dat.p, unit.E)) = porov(5);
    rock.poro(ismember(G_dat.p, unit.F)) = porov(6);
    rock.poro(ismember(G_dat.p, unit.G)) = porov(7);
    
    % ---------------------------------------------------------------------
    % Permeability
    % ---------------------------------------------------------------------
    % Fit a Kozeny-Carman-type equation to data from:
    %   - Beard & Weyl, AAPG (1973)
    %   - Trevisan et al., IJGGC (2014) 
    % Data:
    % diameter_mm = [0.2135, 0.23, 0.45, 0.6, 0.855, 1.4]; 
    % poro_frac = [0.339, 0.379, 0.354, 0.384, 0.38, 0.396];
    % perm_d = [7, 20.17, 115.51, 151, 302, 1570.5];
    
    % Model result from cftool with this custom model, r2 = 0.971
    perm_model = @(x,y) 1.225e+4*x.^2.*y.^3;
    permv = perm_model(d_mm(2:end), porov(2:end));
    rock.perm = nan(G.cells.num, 1);    % isotropic
    if ~removeS
        rock.perm(ismember(G_dat.p, unit.S)) = 0.1*(nano*darcy);
    end
    rock.perm(ismember(G_dat.p, unit.ESF)) = permv(1)*darcy;
    rock.perm(ismember(G_dat.p, unit.ESFsup)) = permv(1)*darcy;
    rock.perm(ismember(G_dat.p, unit.C)) = permv(2)*darcy;
    rock.perm(ismember(G_dat.p, unit.D)) = permv(3)*darcy;
    rock.perm(ismember(G_dat.p, unit.E)) = permv(4)*darcy;
    rock.perm(ismember(G_dat.p, unit.F)) = permv(5)*darcy;
    rock.perm(ismember(G_dat.p, unit.G)) = permv(6)*darcy;
    
elseif model_case == 2 || model_case == 3    % Measurements reported in description.pdf
    rock.poro = nan(G.cells.num, 1);
    rock.perm = nan(G.cells.num, 1);
    
    if ~removeS, rock.poro(ismember(G_dat.p, unit.S)) = .01; end
    rock.poro(ismember(G_dat.p, unit.ESF)) = .435;
    rock.poro(ismember(G_dat.p, unit.ESFsup)) = .435;
    rock.poro(ismember(G_dat.p, unit.C)) = .435;
    rock.poro(ismember(G_dat.p, unit.D)) = .44;
    rock.poro(ismember(G_dat.p, unit.E)) = .45;
    rock.poro(ismember(G_dat.p, unit.F)) = .44;
    rock.poro(ismember(G_dat.p, unit.G)) = .45;
    
    if ~removeS, rock.perm(ismember(G_dat.p, unit.S)) = 0.1*(nano*darcy); end
    rock.perm(ismember(G_dat.p, unit.ESF)) = m(1)*44*darcy;
    rock.perm(ismember(G_dat.p, unit.ESFsup)) = m(2)*44*darcy;
    rock.perm(ismember(G_dat.p, unit.C)) = m(3)*473*darcy;
    rock.perm(ismember(G_dat.p, unit.D)) = m(4)*1110*darcy;
    rock.perm(ismember(G_dat.p, unit.E)) = m(5)*2005*darcy;
    rock.perm(ismember(G_dat.p, unit.F)) = m(6)*4259*darcy;
    rock.perm(ismember(G_dat.p, unit.G)) = m(7)*9580*darcy;
end

if plotPoro == 1
    plts.fig3D(); plotCellData(G, rock.poro, 'edgealpha', 0.2)
    plts.setAxProps(gca), colormap(copper), c = colorbar; 
    if ~removeS, caxis([0.35 .45]), end
    set(get(c,'label'),'string','$\phi$ [-]', 'Interpreter', 'latex', ...
        'fontSize', 14);
    axis equal off; view([90 0]),
end

if plotPerm == 1
    plts. fig3D(); 
    plotToolbar(G, log10(rock.perm/(darcy)), 'edgecolor', [0.5 0.5 0.5], 'edgealpha', 0.4)
    plts.setAxProps(gca), colormap(copper), c = colorbar; 
    if ~removeS, caxis([-2 4]), end
    set(get(c,'label'),'string','$\log_{10}k$ [D]', 'Interpreter', 'latex', ...
        'fontSize', 14);
    axis equal off; view([90 0]), %ylim([0.3 0.6]); zlim([0.42 0.47])
end


end
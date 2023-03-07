function rock = getRock(inj_type, model_case, G, G_dat, unit, m, plts, plotPoro, plotPerm)
%
%
%

%       ESF,  C,   Cf,  E,  Fsup and med, Finf
%d_mm = [0.2, 0.66, 0.66, 1.45, 1.77, 1.77]; % average sand-grain diameters [mm]
%       ESF,  C,   Cf,  E,  Fsup, Finf, Fmed
% d_mm = [0.2, 0.66, 0.66, 1.45, 1.77, 1.77, 1.77];
%d_mm = [0.2, 0.66, 0.66, 1.45, 2.2, 2.2, 2.2];
if model_case == 1    
    % --------------------------------------------------------------------
    % Porosity
    % --------------------------------------------------------------------
    % Data from Beard & Weyl, AAPG (1973) for moderately to well
    % sorted and from Smits et al. (2010).
    rock.poro = nan(G.cells.num, 1);
    if inj_type < 3
        porov = [.37, .38, .38, .39, .39, .39, .39];    % ESF, C, E, F 
        d_mm = [0.1, 0.66, 0.2, 1.45, 3, 3, 2];
        
        rock.poro(ismember(G_dat.p, unit.ESF)) = porov(1);
        rock.poro(ismember(G_dat.p, unit.C)) = porov(2);
        rock.poro(ismember(G_dat.p, unit.Cf)) = porov(3);
        rock.poro(ismember(G_dat.p, unit.E)) = porov(4);
        rock.poro(ismember(G_dat.p, unit.Fsup)) = porov(5);
        rock.poro(ismember(G_dat.p, unit.Finf)) = porov(6);
        rock.poro(ismember(G_dat.p, unit.Fmid)) = porov(7);
    elseif inj_type == 3
        porov = [.37, .38, .38, .38, .39, .39, .39, .39];    % ESF, C, E, F 
        d_mm = [0.1, 0.66, 0.66, 0.66, 1.45, 1.45, 3, 3];
        
        rock.poro(ismember(G_dat.p, unit.ESF)) = porov(1);
        rock.poro(ismember(G_dat.p, unit.Csup)) = porov(2);
        rock.poro(ismember(G_dat.p, unit.Cinf)) = porov(3);
        rock.poro(ismember(G_dat.p, unit.CESF)) = porov(4);
        rock.poro(ismember(G_dat.p, unit.Esup)) = porov(5);
        rock.poro(ismember(G_dat.p, unit.Einf)) = porov(6);
        rock.poro(ismember(G_dat.p, unit.Fsup)) = porov(7);
        rock.poro(ismember(G_dat.p, unit.Finf)) = porov(8);
    end
    
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
    permv = perm_model(d_mm, porov);
    rock.perm = nan(G.cells.num, 1);    % isotropic
    if inj_type < 3
        rock.perm(ismember(G_dat.p, unit.ESF)) = permv(1)*darcy;
        rock.perm(ismember(G_dat.p, unit.C)) = permv(2)*darcy;
        rock.perm(ismember(G_dat.p, unit.Cf)) = permv(3)*darcy;
        rock.perm(ismember(G_dat.p, unit.E)) = permv(4)*darcy;
        rock.perm(ismember(G_dat.p, unit.Fsup)) = permv(5)*darcy;
        rock.perm(ismember(G_dat.p, unit.Finf)) = permv(6)*darcy;
        rock.perm(ismember(G_dat.p, unit.Fmid)) = permv(7)*darcy;
    elseif inj_type == 3
        rock.perm(ismember(G_dat.p, unit.ESF)) = permv(1)*darcy;
        rock.perm(ismember(G_dat.p, unit.Csup)) = permv(2)*darcy;
        rock.perm(ismember(G_dat.p, unit.Cinf)) = permv(3)*darcy;
        rock.perm(ismember(G_dat.p, unit.CESF)) = permv(4)*darcy;
        rock.perm(ismember(G_dat.p, unit.Esup)) = permv(5)*darcy;
        rock.perm(ismember(G_dat.p, unit.Einf)) = permv(6)*darcy;
        rock.perm(ismember(G_dat.p, unit.Fsup)) = permv(7)*darcy;
        rock.perm(ismember(G_dat.p, unit.Finf)) = permv(8)*darcy;
    end

elseif model_case == 2 || model_case == 3          % Measurements reported in description.pdf
    rock.poro = nan(G.cells.num, 1);
    rock.perm = nan(G.cells.num, 1);
    if inj_type < 3
        rock.poro(ismember(G_dat.p, unit.ESF)) = .435;
        rock.poro(ismember(G_dat.p, unit.C)) = .435;
        rock.poro(ismember(G_dat.p, unit.Cf)) = .435;
        rock.poro(ismember(G_dat.p, unit.E)) = .45;
        rock.poro(ismember(G_dat.p, unit.Fsup)) = .44;
        rock.poro(ismember(G_dat.p, unit.Finf)) = .44;
        rock.poro(ismember(G_dat.p, unit.Fmid)) = .44;

        rock.perm(ismember(G_dat.p, unit.ESF)) = m(1)*44*darcy;
        rock.perm(ismember(G_dat.p, unit.C)) = m(2)*473*darcy;
        rock.perm(ismember(G_dat.p, unit.Cf)) = m(3)*473*darcy;
        rock.perm(ismember(G_dat.p, unit.E)) = m(4)*2005*darcy;
        rock.perm(ismember(G_dat.p, unit.Fsup)) = m(5)*4259*darcy;
        rock.perm(ismember(G_dat.p, unit.Finf)) = m(6)*4259*darcy;
        rock.perm(ismember(G_dat.p, unit.Fmid)) = m(7)*4259*darcy;    
    elseif inj_type == 3
        rock.poro(ismember(G_dat.p, unit.ESF)) = .435;
        rock.poro(ismember(G_dat.p, unit.Csup)) = .435;
        rock.poro(ismember(G_dat.p, unit.Cinf)) = .435;
        rock.poro(ismember(G_dat.p, unit.CESF)) = .435;
        rock.poro(ismember(G_dat.p, unit.Esup)) = .45;
        rock.poro(ismember(G_dat.p, unit.Einf)) = .45;
        rock.poro(ismember(G_dat.p, unit.Fsup)) = .44;
        rock.poro(ismember(G_dat.p, unit.Finf)) = .44;

        rock.perm(ismember(G_dat.p, unit.ESF)) = m(1)*44*darcy;
        rock.perm(ismember(G_dat.p, unit.Csup)) = m(2)*473*darcy;
        rock.perm(ismember(G_dat.p, unit.Cinf)) = m(3)*473*darcy;
        rock.perm(ismember(G_dat.p, unit.CESF)) = m(4)*473*darcy;
        rock.perm(ismember(G_dat.p, unit.Esup)) = m(5)*2005*darcy;
        rock.perm(ismember(G_dat.p, unit.Einf)) = m(6)*2005*darcy;
        rock.perm(ismember(G_dat.p, unit.Fsup)) = m(7)*4259*darcy;
        rock.perm(ismember(G_dat.p, unit.Finf)) = m(8)*4259*darcy;
    end
end

if plotPoro == 1
    latx = {'Interpreter', 'latex'};
    fontsz = {'fontSize', 16};
    plts.fig3D(); plotCellData(G, rock.poro, 'edgealpha', 0.2)
    plts.setAxProps(gca), colormap(copper), c = colorbar; 
    set(get(c,'label'),'string','$\phi$ [-]', 'Interpreter', 'latex', ...
        'fontSize', 14);
    axis equal off; view([90 0]),

    f= figure(11);
    h = subplot(1,2,1);
    plot(1:numel(d_mm), porov, 'sk', 'markerSize', 10, ...
         'markerfacecolor', [0.3 0.3 0.3])
    grid on
    xlabel('Sand type', fontsz{:})
    ylabel('$n$ [-]', latx{:}, fontsz{:})
    xticklabels({'ESF', 'C', 'Cf', 'E', 'Finf', 'F'})
    ylim([0.35 0.4])
    h.FontSize = 14;
    
    h = subplot(1,2,2);
    plot(1:numel(d_mm), log10(permv), 'sk', 'markerSize', 10, ...
        'markerfacecolor', [0.3 0.3 0.3])
    grid on
    xlabel('Sand type', fontsz{:})
    ylabel('$\log_{10}(k$ [D])', latx{:}, fontsz{:})
    ylim([0 4])
    xticklabels({'ESF', 'C', 'Cf', 'E', 'Finf', 'F'})
    h.FontSize = 14;
    f.Position = [0, 0, 600, 250];
end

if plotPerm == 1
    plts. fig3D(); plotCellData(G, log10(rock.perm/(darcy)), 'edgealpha', 0.2)
    plts.setAxProps(gca), colormap(copper), c = colorbar; %caxis([3 4])
    set(get(c,'label'),'string','$\log_{10}k$ [D]', 'Interpreter', 'latex', ...
        'fontSize', 14);
    axis equal off; view([90 0]), %ylim([0.3 0.6]); zlim([0.42 0.47])
end


end
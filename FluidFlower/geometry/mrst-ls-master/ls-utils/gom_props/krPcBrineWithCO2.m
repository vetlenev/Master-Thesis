function krPcBrineWithCO2(T, P, opt, figs, directory)
%
% Lluis Salo, October 2019
% lsalo@mit.edu
%
%

%% Check and organize inputs
% ------------------------- Temperature -----------------------------------
T_unit = T{1}; T_val  = T{2};
if strcmp(T_unit, 'K'),     T_K = T_val;
elseif strcmp(T_unit, 'C'), T_K = T_val + 273.15;
elseif strcmp(T_unit, 'F'), T_K = (T_val + 459.67)*(5/9);
else, error("Temperature units must be 'K', 'C' or 'F'.")
end
T_C   = T_K - 273.15;                                                                                                    

% ---------------------------- Pressure -----------------------------------
p_unit = P{1}; p_val  = P{2};
assert(max(strcmp(p_unit, {'Pa','MPa','bar','psi'})))
if strcmp(p_unit, 'Pa'),        P = p_val*10^(-5);
elseif strcmp(p_unit, 'MPa'),   P = p_val*10;
elseif strcmp(p_unit, 'psi'),   P = p_val*0.06894757;
elseif ~strcmp(p_unit, 'bar'), error("Pressure units not supported")
end
%P_kPa = P*100;


%% kR (Clay-rich and sand-rich)
switch opt.kr
    case 'water-co2-ExxonComparison'
        Swr = 0.56;
        S_gr = 0;
        m    = 3;
        n    = 3;
        krw0 = 1;
        krg0 = 0.33;
        S    = fliplr(linspace(Swr, 1-S_gr,31));
        Sg   = 1 - S; 
        S    = [S(1) 0.999 S(2:end)];  Sg = [Sg(1) 0.001 Sg(2:end)];   
        kr_w = krw0.*(1 - (Sg./(1-Swr))).^n;
        kr_n = krg0.*(Sg./(1-Swr)).^m;
        
    case 'Miocene-base'
        nreg  = 2;
        %poro  = [0.3 0.15];
        %perm  = [150, 0.01];
        Swr  = [0.3092; 0.43];
        f     = [0.7; 0.4];
        S_gr  = zeros(nreg,1);
        n     = 2.6;
        S     = fliplr([linspace(Swr(1),1-S_gr(1),32); linspace(Swr(2),1-S_gr(2),32)]);
        Sg    = 1 - S;  Sg(:, 2) = 0.001;   S(:, 2) = 0.999;
        nS    = (S-Swr) ./ (1-Swr-S_gr);
        kr_w  = nS.^n;
        kr_n  = f.*(1-nS).^n;
        
        Sgt  = [0.18; 0.3];
        Sgi   = [0 linspace(Sgt(1), 1-Swr(1), 31); 0 linspace(Sgt(2), 1-Swr(2), 31)];
        Swi   = 1 - Sgi;
        nSi   = (Swi-Swr) ./ (1-Swr-Sgt); nSi(nSi<0)=0;
        kr_wi = [interp1(S(1,:), kr_w(1,:), Swi(1,:)); interp1(S(2,:), kr_w(2,:), Swi(2,:))];
        kr_wi(:,end) = 0;
        kr_ni = f.*(1-nSi).^n; kr_ni(:,1:2) = 0;
        
    case 'Miocene-case2' % Same as above but with 3rd region for fault
                         % Now Land's model is used to calculate bounding
                         % imbition curve
        nreg  = 3;
        %poro  = [0.3 0.15];
        %perm  = [150, 0.01];
        Swr  = [0.3092; 0.43; 0.3696];
        f     = [0.9; 0.6; 0.8]; % krg at Sgmax
        S_gr  = zeros(nreg,1);
        e      = 4;                  % typical for water in sandstones, using Sw*
        lambda = 2/(e-3);           % Land (1968)
        ng     = 2/lambda + 1;       % Brooks & Corey (1964), Eq. 15
        S     = fliplr([linspace(Swr(1),1-S_gr(1),32); 
                        linspace(Swr(2),1-S_gr(2),32);
                        linspace(Swr(3),1-S_gr(3),32);]);
        Sg    = 1 - S;  Sg(:, 2) = 0.001;   S(:, 2) = 0.999;
        nS    = (S-Swr) ./ (1-Swr-S_gr);
        kr_w  = nS.^e;
        kr_n  = f.*(1-nS).^ng;              % Simplified
        kr_n2 = f.*((1-nS).^2).*(1-nS.^ng); % Actual Brooks & Corey for drainage
        
        % Imbibition
        Sgt  = [0.18; 0.3; 0.24];
        Sgi   = [linspace(Sgt(1), 1-Swr(1), 32); linspace(Sgt(2), 1-Swr(2), 32); ...
                 linspace(Sgt(3), 1-Swr(3), 32)];
        Swi   = 1 - Sgi;
        
        % Land-based bounding imbibition curve (Land, 1968)
        % Gas 
        nsgt_max = Sgt ./ (1-Swr);
        C = 1./nsgt_max - 1;                           % Land (1968) Eq. 2
        nsgib = (Sgi)./(1-Swr);  
        for n=1:nreg                                   % Land (1968) Eq. 4
            nsgibF(n,:) = 0.5*((nsgib(n,:)-nsgib(n,1)) + ...
                               sqrt((nsgib(n,:)-nsgib(n,1)).^2 + ...
                                    (4/C(n))*(nsgib(n,:)-nsgib(n,1))));  
        end
        %e = 4;                                        % Land (1968) Fig. 2
        krgib = nsgibF.^2.*(1-(1-nsgibF).^(e-2));
        kr_ni = f.*krgib;                              % our krg is not 1 at sgi;
        
        % Water (no hysteresis, just interpolate at new saturations)
        nSi   = (Swi-Swr) ./ (1-Swr-Sgt); nSi(nSi<0)=0;
        kr_wi = [interp1(S(1,:), kr_w(1,:), Swi(1,:)); interp1(S(2,:), kr_w(2,:), Swi(2,:)); ...
                 interp1(S(3,:), kr_w(3,:), Swi(3,:))];
        kr_wi(:,end) = 0;
        % Old gas imbibition curve
        % kr_ni = f.*(1-nSi).^n; kr_ni(:,1:2) = 0;
    otherwise
        error('kr case not supported')
end


%% Pc (Clay-rich)
switch opt.Pc
    case 'Miocene-base'
        % Data for MICP curves from Lu et al. (2017), which are approximately
        % averaged. Data for contact angles and interfacial tensions are 
        % for T = 70C, P=20MPa.
        micp_raw    = [0, 1000, 1250, 1500, 1800, 2100, 2550, 3000, 3650, ...
                       4300, 5400, 6500, 8000, 9500, 11250, 13000, ...
                       16000, 19000];
        sg0         = 0:0.05:0.8; sg0 = [sg0(1) 0.0001 sg0(2:end)];
        micp        = interp1(sg0, micp_raw, Sg(2,:));
        theta       = [70, 140];
        ift         =  [25, 485];
        psi_to_bar  = 0.06894757;
        Pc_co2br(2,:) = psi_to_bar*micp*abs( cosd(theta(1))*ift(1)/(cosd(theta(2))*ift(2)) );
        Pc_co2br(1,:) = Pc_co2br(2,:).*0.01; % JUST FOR TEST, this sat region should be left 0
        
        sgii           = [0, 0.3, 0.31, 0.35, 0.4, 0.45, 0.5, 0.55, 0.57];
        pcii           = [0,   0,  0.9,  1.2, 1.6,  2.5, 5.0,  9.0, Pc_co2br(2,end)];
        pcii           = interp1(sgii, pcii, Sgi(2,:)); pcii(end) = Pc_co2br(2,end);
        Pci_co2br(2,:) = pcii;
        Pci_co2br(1,:) = Pci_co2br(2,:).*0.01; % JUST FOR TEST
    case 'Miocene-definitive' % Lu et al., 2017 (Green curve)        
        micp          = [700, 850, 1100, 1400, 2000, 2400, 3000, 3500, 4500, ...
                         5700, 7000, 8500, 9600, 10500, 12000, 14000 18000];
        sg_micp       = [0.001 0.05:.05:0.8];
        theta         = [70, 140];
        ift           = [25, 485];
        psi_to_bar    = 0.06894757;  
        Pc_co2br_raw  = psi_to_bar*micp*abs( cosd(theta(1))*ift(1)/(cosd(theta(2))*ift(2)) );  % [bar] 
        
        % Brooks and Corey (1964) model fit
        %pc_e          = Pc_co2br_raw(1);                                    % entry pc
        %Sw            = 1 - sg_micp;
        %nSat          = (Sw-S_wr(3)) ./ (1-S_wr(3)-S_gr(3));
        %idNeg         = nSat<0;
        %nSat(idNeg)   = [];
        %pcToFit       = Pc_co2br_raw;   pcToFit(idNeg) = [];
        %ft            = fittype('pc_e*(x)^(-1/lambda)', 'independent', 'x', ...
        %                        'coefficients', 'lambda', 'problem', 'pc_e');
        %coreyModel    = fit(nSat', pcToFit', ft, 'problem', pc_e);  
        
        % van Genuchten model fit with cftool
        % pg = 4.675;    m = 0.6992;    R^2 = 0.91
        %vGmodel = @(Sg, S_wr) real(4.675*((((1-Sg)-S_wr) ./ (1-S_wr))^(-1/0.6992) - 1)^(1-0.6992));
        
        % Since the BC model (power-law) doesn't fit that well, and the vG
        % model poses convergence problemes and is discouraged, we use
        % linear interpolation on those datapoints.
        Pc_co2br_raw = [0 Pc_co2br_raw];
        sg_micp = [0 sg_micp];
        Pc_co2br(2,:)  = interp1(sg_micp, Pc_co2br_raw, Sg(2,:));
        Pc_co2br(3,:)  = interp1(sg_micp, Pc_co2br_raw, Sg(3,:));
        
        sgii           = [0, 0.3, 0.31, 0.35, 0.4, 0.45, 0.5, 0.55, 0.57];
        pcii           = [0,   0,  0.9,  1.2, 1.6,  2.5, 5.0,  9.0, Pc_co2br(2,end)];
        Pci_co2br(2,:) = interp1(sgii, pcii, Sgi(2,:)); Pci_co2br(2, end) = Pc_co2br(2,end);
        sgii = [0, 0.24, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.57, Sgi(3, end)]; 
        pcii = [0,    0,  0.6, 0.7, 1,    1.3, 1.9,  3, 5,  9,  Pc_co2br(3,end)];
        Pci_co2br(3,:) = interp1(sgii, pcii, Sgi(3,:)); Pci_co2br(3, end) = Pc_co2br(3,end);
        Pci_co2br(1,:) = 0; 
        
    otherwise
        %error('Pc case not supported')
end

%% Plots
if figs == 1
    figure(1)
    if nreg == 2
        subplot(1,2,1)
        plot(S(1,:), kr_w(1,:), '-b'); grid on
        hold on; plot(S(1,:), kr_n(1,:), '-k'); 
        plot(Swi(1,:), kr_ni(1,:), '--', 'color', [0.5 0.5 0.5]);
        %plot(Swi(1,:), kr_wi(1,:), '--c');
        xlim([0 1]); hold off
        xlabel('S_w'), ylabel('k_r Reservoir');
        %legend('$k_{r,w}$', '$k_{r,n}^d$', '$k_{r,n}^i$', 'Interpreter', 'latex', 'fontsize', 12)
        set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1])
        set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1])
        subplot(1,2,2)
        plot(S(2,:), kr_w(2,:), '-b'); grid on
        hold on; plot(S(2,:), kr_n(2,:), '-k'); 
        plot((1-Sgi(2,:)), kr_ni(2,:), '--', 'color', [0.5 0.5 0.5]);
        %plot(Swi(2,:), kr_wi(2,:), '--c');
        xlim([0 1]); hold off
        xlabel('S_w'), ylabel('k_r Seal');
        legend('$k_{r,w}$', '$k_{r,n}^d$', '$k_{r,n}^i$', 'Interpreter', 'latex', 'fontsize', 12)
        set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1])
        set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1])
    elseif nreg == 3
        latx = {'interpreter','latex'};
        h = figure(1);
        tiledlayout(1,3,'TileSpacing','Compact');
        nexttile
        p1 = plot(S(1,:), kr_w(1,:), '-b', 'linewidth', 1, ...
                  'displayname', '$k_{\mathrm{r,aq}}$'); 
        grid on
        hold on; 
        p2 = plot(S(1,:), kr_n2(1,:), '-r', 'linewidth', 1, ...
                  'displayname', '$k_{\mathrm{rg}}$'); 
        %plot(S(1,:), kr_n(1,:), '.r'); 
        p3 = plot(Swi(1,:), kr_ni(1,:), '--r', 'linewidth', 1, ...
                  'displayname', '$k_{\mathrm{rg}}^\mathrm{i}$');
        xlim([0 1]); hold off
        xlabel('$S_\mathrm{w}$ [-]','fontsize',12, latx{:})
        ylabel('$k_\mathrm{r}$ [-]','fontsize',12, latx{:})
        legend([p1 p2 p3], latx{:},'fontsize',11)
        set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1])
        set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1])
        title('Storage reservoir (SR)', 'fontweight', 'normal')
        
        nexttile
        colororder({'k','k'})
        yyaxis left
        plot(S(2,:), kr_w(2,:), '-b', 'linewidth', 1); grid on
        hold on; 
        plot(S(2,:), kr_n2(2,:), '-r', 'linewidth', 1); 
        %plot(S(2,:), kr_n(2,:), '.r'); 
        %plot(Swi(2,:), kr_ni(2,:), '--r');
        set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1])
        ylabel('$k_\mathrm{r}$ [-]','fontsize',12, latx{:})
        yyaxis right
        p1 = plot(1-Sg(2,:), Pc_co2br(2,:), '-k', 'linewidth', 1, ...
                  'displayName', '$P_\mathrm{c}$');
        ylabel('$P_\mathrm{c}$ [bar]', latx{:}, 'fontsize', 12),
        xlim([0 1]); hold off
        set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1])
        set(gca,'Ytick', 0:2:14)
        xlabel('$S_\mathrm{w}$ [-]','fontsize',12, latx{:})
        title('Top seal (TS)', 'fontweight', 'normal')
        legend(p1, latx{:}, 'fontsize', 11)
        
        nexttile
        plot(S(3,:), kr_w(3,:), '-b', 'linewidth', 1); grid on
        hold on; plot(S(3,:), kr_n2(3,:), '-r', 'linewidth', 1);
        %plot(S(3,:), kr_n(3,:), '.r');
        %plot(Swi(3,:), kr_ni(3,:), '--r');
        xlim([0 1]); hold off
        set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1])
        set(gca,'Ytick',[0 0.2 0.4 0.6 0.8 1])
        ylabel('$k_\mathrm{r}$ [-]','fontsize',12, latx{:})
        xlabel('$S_\mathrm{w}$ [-]','fontsize',12, latx{:})
        title('Fault', 'fontweight', 'normal')
        set(h, 'position', [50, 50, 800, 250])
    end
    
    h = figure(2);
    plot(1-Sg(2,:), Pc_co2br(2,:), '-k', 'linewidth', 1);
    grid on, xlim([0 1]);
    xlabel('$S_w$ [-]', latx{:}, 'fontsize', 12) 
    ylabel('$P_\mathrm{c}$ [bar]', latx{:}, 'fontsize', 12),
    title('Top seal (TS)', 'fontweight', 'normal')
    set(h, 'position', [50, 50, 400, 350])
    
%     figure(2)
%     subplot(1,2,1)
%     plot(1-Sg(2,:), Pc_co2br(2,:), '-xk', 'linewidth', 1.5); hold on
%     plot(1-Sgi(2,:), Pci_co2br(2,:), 'o-', 'color', [0.7 0.7 0.7], 'linewidth', 1.5);
%     grid on, xlim([0 1]);
%     xlabel('S_w'), ylabel('P_c [bar]'), legend('Drainage', 'Imbibition')
%     title('Amph b')
%     subplot(1,2,2)
%     plot(1-Sg(3,:), Pc_co2br(3,:), '-xk', 'linewidth', 1.5); hold on
%     plot(1-Sgi(3,:), Pci_co2br(3,:), 'o-', 'color', [0.7 0.7 0.7], 'linewidth', 1.5);
%     grid on, xlim([0 1]);
%     xlabel('S_w'), ylabel('P_c [bar]'), legend('Drainage', 'Imbibition')
%     title('Fault (max Vcl)')
    
end

%% Write table
% Write to file
if nargin > 4  
    fileName  = opt.kr;
    finImb    = 'SGOF_imb';
    if nreg == 2
        satg = [Sg(1,:) NaN Sg(2,:)]';
        krn =  [kr_n(1,:) NaN kr_n(2,:)]';
        krw =  [kr_w(1,:) NaN kr_w(2,:)]';
        Pc   = [Pc_co2br(1,:) NaN Pc_co2br(2,:)]';
        
        sgi   = [Sgi(1,:) NaN Sgi(2,:)]';
        krwi  = [kr_wi(1,:) NaN kr_wi(2,:)]';
        krni  = [kr_ni(1,:) NaN kr_ni(2,:)]';
        Pci   = [Pci_co2br(1,:) NaN Pci_co2br(2,:)]'; 
    elseif nreg == 3
        satg = [Sg(1,:) NaN Sg(2,:) NaN Sg(3,:)]';
        krn = [kr_n2(1,:) NaN kr_n2(2,:) NaN kr_n2(3,:)]';
        krw = [kr_w(1,:) NaN kr_w(2,:) NaN kr_w(3,:)]';
        Pc   = [Pc_co2br(1,:) NaN Pc_co2br(2,:) NaN Pc_co2br(3,:)]';
        
        sgi   = [Sgi(1,:) NaN Sgi(2,:) NaN Sgi(3,:)]';
        krwi  = [kr_wi(1,:) NaN kr_wi(2,:) NaN kr_wi(3,:)]';
        krni  = [kr_ni(1,:) NaN kr_ni(2,:) NaN kr_ni(3,:)]';
        Pci   = [Pci_co2br(1,:) NaN Pci_co2br(2,:) NaN Pci_co2br(3,:)]';
    end
    tab = table(satg, krn, krw, Pc);
    tab.Properties.VariableNames = {'SGAS' 'KRG' 'KRWG' 'PCOG'};
    writetable(tab, [directory fileName '.txt'], 'Delimiter', 'tab');
    
    tabi = table(sgi, krni, krwi, Pci);
    tabi.Properties.VariableNames = {'SGASi' 'KRGi' 'KROGi' 'PCOGi'};
    writetable(tabi, [directory finImb '.txt'], 'Delimiter', 'tab');
end

end
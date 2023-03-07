%% Pce as a fcn of perm-poro
% Since Leverett scaling on the basis of one reference curve does not 
% seem to work well for clays, we model Pce as a fcn of porosity and
% permeability using data in Treviño & Meckel, GoM Atlas (2017); Ewy,
% GEE 2018; Tokunaga, WRR 2020
% From Pore radius Fig 11 in Ewy et al MPG (2020) we get Pce Hg-Air between
% 92 bar (1340 psi) (shale A, depth 1000) and 185 bar (2690 psi) 
% (shale G-H, depth 2000-3000), which is within range of GoM samples,
% although perm is in the nD range here.
% Pce = 2*ift*cos(theta)/max_pore_radius

% Data     [[GoM Atlas] [Ewy] [Tokunaga2017(4) Al-Bazali 2005(3)]]
porosity = [[6.48 7.81 5.65 8.36 7.04 3.15]*0.01 ...
            [0.4 0.22 0.15]./[1.4 1.22 1.15] [8.1 7 7 10 nan nan nan]*0.01];
k_md = [[1.64 1.72 0.92 2.05 2.08 0.145]*1e-3 ...
        [7.5 0.3 0.03]*1e-21/(milli*darcy) ...
        [0.04 5 7 0.01 0.65 0.045 0.296]*1e-20/(milli*darcy)];
watAir_to_HgAir = abs(cosd(140)*485/(cosd(5)*72.8)); % myNotes for references
psi_to_bar = 0.06894757;
pce_bar = [[1303, 1539, 1200, 520, 613, 2535]*psi_to_bar ...
           [1.8 3.5 3.5]*10*watAir_to_HgAir ...
           [100 9 9 9 43.4 65.5 48.3]*watAir_to_HgAir]; % All converted to Hg-Air
weights=[1 1 1 1 1 1 0.5 0.5 0.5 0.1 0.25 0.25 0.1 0.25 0.25 0.25];

% From cftool with k_md using a linear model (not weighted)
% This works in the range of the data above, i.e. log10(k_md) \in [-7.5, -2.5]
% R2 = 0.28 though...
pce_model_hg = @(perm_md) -40.85*log10(perm_md) -41.19;     % bar
pce_busch13 = @(k_m2) -0.11*log10(k_m2) - 1.85;               % log MPa

% to convert to co2-brine
pce_model_co2brine = @(pce_hg, theta) pce_hg*abs(cosd(theta)*25/(cosd(140)*485));

%% GoM-consistent data only
% Tertiary sediments + clay smear Sperrevik (3 samples). We take Pd 10% 
% where curve is available or where reported (GoM Atlas, AH, H), 
% displacement pressure (B-A) and threshold pressure (Sperrevik 2002).
% Since the correlation is pretty poor for broad datasets (see Fig. 10 and
% 11 only mudrocks in Busch & AH MPG 2013) we picked a smaller dataset of
% geologically more consistent data. The correlation is still weak but
% better, and this incrases confidence in our predictions (still consider
% +- RMSE as possible). The use of displacement pressure is better over
% entry pressure because in continuum reservoir models (coarse) the first 
% nonzero Pc value allows flow in a large (tens of m) section.
% *******************
% TBD: Update GoM and Heath et al to Pd,10, then adjust final model and move on.
% *******************
k_md2 = [[1.64 1.72 0.92 2.05 2.08 0.145]*1e-3 ...               % Mudrock Miocene GoM Atlas
         [36.5 1.6]*1e-3 ...                                     % D-A 2002 (Miocene Deep Water GoM) (excluded superseal!)
         [3.4 mean([34.9 17.1]) 9.9]*nano*darcy/(milli*darcy) ...% H et al 2004 (M2-M4; Tertiary mudstone NOR Shelf)
         [0.05 0.05 0.9 0.9]*1e-20/(milli*darcy) ...             % Heath et al. Tuscaloosa Late Cret GoM
         [0.09 0.07 0.09 0.12] ...                               % B-A 1995 (Stegall 1A, Tertiary growth faults Texas Gulf Coast)
         [5 3 1]*1e-5];                                          % Sperrevik 02 clay smears
pce_bar2 = [[1850, 2200, 1800, 1100, 1050, 3400]*psi_to_bar ...
            [465 4845]*psi_to_bar ...
            [23 2.6 16.2]*10 ...
            [66 29 56 66]*10 ...
            [547 554 544 410]*psi_to_bar ...
            [2000 3500 3500]*psi_to_bar];
pcd_fit  = fittype('a*x + b', 'independent', 'x');
[pcd_model, g] = fit(log10(k_md2)', log10(pce_bar2)', pcd_fit);
[rho, pval] = corr(log10(k_md2)', log10(pce_bar2)');

% Busch 2013 Fig 10 correlation
pcb_busch13 = @(k_m2) -0.25*log10(k_m2) -4.337;      % log MPa
watHe_to_HgAir = abs(cosd(140)*485/(cosd(0)*70));    % He, 25C, 0-50 MPa


%% Clay content (mineralogy) vs Pc (only GoM)
% Data from Lu et al. (2017) and Dawson & Almon (2002), Dawson & Almon
% (2006), Day-Stirrat (2012)
% Turns out NOT to be better than correlating for permeability.
cc = [[31.9 38.1 35.8 30.6 32.8 27.3] ... % GoM Atlas
      [52 59 54 59 64 67 72 69 72] ...    % Miocene Deepwater GoM
      [76 69 61 57 49 64] ...             % Miocene Deepwater GoM
      [53.7 52.1 34.4 41.5]];             % Ursa basin (< 400 mbsf) Quaternary GoM
pce_bar3 = [[1303, 1539, 1200, 520, 613, 2535] ...
            [465 2530 4845 4645 5665 6260 6775 8155 18195] ...
            [8394 7448 4511 3177 1360 7653] ...
            [570 365 1148 1240]]*psi_to_bar;
            

%% Plots

% To calculate Pc(sg, log10_perm_md)
% pcv(sg) = ( pce_model_co2brine(pce_model_hg(log10_perm_md),
%                                theta_val)/refPc(sg_entry) )*refPc(sg)

latx = {'interpreter', 'latex'};
f = figure(141);
hold on
p1 = plot(log10(k_md), pce_bar, 'sk', 'DisplayName', 'data');
vlogk = -8:.1:-2;
p2 = plot(vlogk, pce_model_hg(10.^vlogk), '-k', 'linewidth', 1, ...
          'DisplayName', '$P(x)$ fit');
p3 = plot(vlogk, 10*(10.^pce_busch13(10.^vlogk*(milli*darcy))), '-b', ...
          'linewidth', 1, 'DisplayName', 'B-AH 13 fit');
text(-4, 400, 'R$^2 = 0.28$', latx{:}, 'fontSize', 10)
hold off
grid on
xlabel('$\log_{10}(k$ [mD])', latx{:}, 'fontsize', 12)
ylabel('$P_c^\mathrm{e}$ [bar]', latx{:}, 'fontsize', 12)
xlim([-8 -2]), xticks(-8:1:-2)
ylim([0 550]), yticks(0:100:500)
h=legend([p1 p2 p3], 'location', 'northeast', latx{:}, 'fontsize', 10);
set(h.BoxFace, 'ColorType','truecoloralpha', ...
            'ColorData', uint8(255*[1;1;1;.5])); 
f.Position = [100 100 350 300];

f = figure(142);
hold on
p1 = plot(log10(k_md2(1:6)), pce_bar2(1:6), 'sk', 'DisplayName', '$P_{d,10}$ Lu 2017 (Miocene TX Shelf), $z \sim 3.2$km', 'markerFaceColor', 'k');
p2 = plot(log10(k_md2(7:8)), pce_bar2(7:8), 'sk', 'DisplayName', '$P_{d,10}$ Dawson 2002 (Miocene Deep Water GoM)', 'markerFaceColor', [0.6 0.6 0.6]);
p3 = plot(log10(k_md2(9:11)), pce_bar2(9:11), 'sk', 'DisplayName', '$P_{d,10}$ Hildenbrand 2004 (Tertiary NOR Shelf), $z \sim 1.5$km', 'markerFaceColor', 'r');
p4 = plot(log10(k_md2(12:15)), pce_bar2(12:15), 'sk', 'DisplayName', ' $P_b$ Heath 2011 (Late Cretaceous Gulf Coast) $z \sim 2.5$km', 'markerFaceColor', 'm');
p5 = plot(log10(k_md2(16:19)), pce_bar2(16:19), 'dk', 'DisplayName', ' $P_d$ Berg 1995 (Shear zone, Eocene Gulf Coast) $z \sim 4.2$km', 'markerFaceColor', 'b');
p6 = plot(log10(k_md2(20:end)), pce_bar2(20:end), 'dk', 'DisplayName', '$P_t$ Sperrevik 2002 (Clay smear, N Sea / NOR Shelf)', 'markerFaceColor', 'c');
vlogk = -7:.1:0;
p7 = plot(vlogk, 10.^(pcd_model(vlogk)), '-k', 'linewidth', 1, ...
          'DisplayName', '$-0.1992x +1.407$, R$^2 = 0.59$, $p = 3\times10^{-5}$');
p8 = plot(vlogk, 10.^(pcd_model(vlogk)-g.rmse), '--k', 'linewidth', 0.5, ...
          'DisplayName', '$\pm$ RMSE');
p9 = plot(vlogk, 10.^(pcd_model(vlogk)+g.rmse), '--k', 'linewidth', 0.5, ...
          'DisplayName', '$P(x)$ fit');
pcb = 10.^(pcb_busch13((10.^vlogk)*milli*darcy))*watHe_to_HgAir*10; % log bar Hg-Air
p10 =  plot(vlogk, pcb, '-', 'color', [0.7 0.7 0.7], 'linewidth', 1, ...
          'DisplayName', 'Busch 2013 (Eq. 27), H$_2$O-He: $\theta=0, \gamma = 70$mN/m');
hold off
grid on
xlabel('$\log_{10}(k$ [mD])', latx{:}, 'fontsize', 12)
ylabel('Hg-Air $P_c$ [bar]', latx{:}, 'fontsize', 12)
xlim([-7 -0]), xticks(-8:1:0)
ylim([10 10^4]), yticks([10 100 1000 10^4])
h=legend([p1 p2 p3 p4 p5 p6 p7 p8 p10], 'location', 'northeast', latx{:}, 'fontsize', 8);
set(h.BoxFace, 'ColorType','truecoloralpha', ...
            'ColorData', uint8(255*[1;1;1;.5]));
set(gca, 'YScale', 'log')
f.Position = [100 100 400 375];
%% Import relPerm data, visualize, compare and pick rock curves for upscaling
%
% Data from ExxonMobil (Lisa S Lun, Oct 2021)
%
close all
mrstModule add ad-props deckformat ad-core ls-proj ls-utils

% Read
%pth = 'C:\Users\lsalo\data_from_Exxon\';
pth = 'C:\Users\Lluis\PhD_MIT\3.research\2020_2021_ExxonMobil\data_from_Exxon\';
%pth = '/media/lsalo/DADES/PhD_MIT/3.research/2020_2021_ExxonMobil/data_from_Exxon/';
%pth = '/Users/lluis/Dropbox (MIT)/PhD_MIT/3.research/2020_2021_ExxonMobil/data_from_Exxon/';
fn = 'WestDelta_30_41_73.xlsx';
data = readtable([pth, fn], 'Sheet', 2);

% Visualize
idA2 = 3:26;
pc_gw_A2 = cell(numel(idA2), 1);
psi_to_pa = 6894.76;
figure(1)
subplot(1,2,1)
hold on
to_res = 25*cosd(30)/72; % [\sigma*cos(\theta)]res / [\sigma*cos(\theta)]lab . 
                         % GoM atlas, for CO2-brine res; Data from Lisa for lab.
for n=idA2
   pc_gw_A2{n-2} = table2array(readtable([pth, fn], 'Sheet', n)); 
   pc_gw_A2{n-2}(:, 2) = to_res*psi_to_pa.*pc_gw_A2{n-2}(:,2);
   
   if n==3
   l1 = plot(pc_gw_A2{n-2}(:, 1), pc_gw_A2{n-2}(:, 2)/1e5, '.-', 'markersize', 4, ...
        'color', [0.7 0.7 0.7], 'DisplayName','1A.2 data');
   else
       plot(pc_gw_A2{n-2}(:, 1), pc_gw_A2{n-2}(:, 2)/1e5, '.-', 'markersize', 4, ...
            'color', [0.7 0.7 0.7])
   end
end
latx  = {'Interpreter','latex'};
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$P_\mathrm{c}$ [bar]', latx{:}, 'fontSize', 14)
title('1A.2: g-w (drainage) PPM $\vert$ z = [1471, 2797]m', latx{:})
%set(gca, 'yscale', 'log')
ylim([0.001 0.7]), yticks([0.01 0.1 0.2 0.3 0.4 0.5 0.6 0.7])
xticks(0:.2:1)
grid on

% Get min, max perm and associated poro
kg.vals = table2array(data(1:numel(idA2), 51))*milli*darcy;
[kg.minv, idmin] = min(kg.vals);
[kg.maxv, idmax] = max(kg.vals);

poro.vals = table2array(data(1:numel(idA2), 48));
poro.minv = min(poro.vals(idmin));
poro.maxv = max(poro.vals(idmax));

% Compare to base-case Pc
fnf = 'ls-proj/gom/adbo_2.5D/eclipse_data_files/krHyst_fPcFault_co2brine.DATA';
deck = convertDeckUnits(readEclipseDeck(fnf));
fluid = initDeckADIFluid(deck);
id_pc_fault = 3;
refPc.val  = fluid.pcOG{id_pc_fault};
refPc.poro = 0.07;                            % Table 3.3, sample 5 in Treviño & Meckel, 2017 (Pc curve is the green one in Fig. 3.4a)
refPc.perm = 0.00208*milli*darcy;
sG = 1-(0.3696:.0025:1);
pcBase.maxv = refPc.val(sG).*sqrt((refPc.perm*poro.minv)./(refPc.poro*kg.minv));
pcBase.minv = refPc.val(sG).*sqrt((refPc.perm*poro.maxv)./(refPc.poro*kg.maxv));

l2 = plot(1-sG, pcBase.minv/1e5, '--r', 'linewidth', 1, ...
          'DisplayName','Scaled from clay-rich core');
plot(1-sG, pcBase.maxv/1e5, '--r', 'linewidth', 1)

% Averaged Pc curve
idrem = 5;      % Curve too different from others
swr = table2array(data(1:numel(idA2), 91));
swr(idrem) = [];
%swr_avg_sand = round(mean(swr), 2);
swr_avg_sand = 0.2;
sg = 0:0.005:(1-swr_avg_sand);
pc_gw_A2_interp = zeros(numel(sg), numel(idA2)+1);
pc_gw_A2_interp(:,1) = 1-sg;
for n=1:numel(pc_gw_A2)
    pc_gw_A2_interp(:,n+1) = interp1(pc_gw_A2{n}(:,1), pc_gw_A2{n}(:,2), 1-sg);
end
pc_avg = zeros(numel(sg), 2);
pc_avg(:, 1) = 1-sg;
ids = 1:numel(idA2);
ids(idrem+1) = [];
for j = 1:numel(sg)
    pc_avg(j, 2) = nanmean(pc_gw_A2_interp(j, ids));
    if sg(j) > 0.6 && pc_avg(j, 2) < pc_avg(j-1, 2)
      pc_avg(j, 2) = 1.02*pc_avg(j-1, 2); 
    end
end
%pc_avg(:, 2) = smooth(pc_avg(:, 1), pc_avg(:, 2), 'rloess');
l3 = plot(pc_avg(:,1), pc_avg(:,2)/1e5, '--', 'linewidth', 1, ...
          'DisplayName', 'mean(1A.2)', 'color', [0.3 0.3 0.3]);
l4 = plot(pc_gw_A2{15}(:,1), pc_gw_A2{15}(:,2)/1e5, '-k', 'linewidth', 2, ...
          'DisplayName', 'A.2.15');


% Fit Brooks-Corey (1964) model
% https://petrowiki.spe.org/Capillary_pressure_models
pc_t          = pc_avg(10,2);                                    % entry pc
Sw            = 1 - sg;
nSat          = (Sw-swr_avg_sand) ./ (1-swr_avg_sand-sg(1));
bc_fit        = fittype('pc_t*(x)^(-1/lambda)', 'independent', 'x', ...
                        'coefficients', 'lambda', 'problem', 'pc_t');
bc_model    = fit(nSat(1:end-1)', pc_avg(1:end-1,2), bc_fit, 'problem', pc_t);
bc_vals = bc_model(nSat(1:end-1));
l5 = plot([1; pc_avg(1:end-1,1)], [1e-5; bc_vals/1e5], '-r', 'linewidth', 2, ...
          'DisplayName', 'BC fit');
 
legend([l1, l3, l2, l4, l5], 'location', 'southwest')

% Plot final Pc curves for sand (this dataset) and for clay material (GoM 
% atlas)
% Final Pc curves in simulation model will be scaled, using the avg poro
% and perm for cores used in this experiment 1A.2 (sands) and the poro and perm
% for the Miocene core (clays).
% Lu et al. (2017) Green curve
micp          = [700, 850, 1100, 1400, 2000, 2400, 3000, 3500, 4500, ...
                 5700, 7000, 8500, 9600, 11000, 13000, 15000 21000];
sg_micp       = [0.001 0.05:.05:0.8];
theta         = [30, 140];  %30
ift           = [25, 485];
fc            = 0.06894757;   % psi to bar
Pc_co2br_raw  = fc*micp*abs( cosd(theta(1))*ift(1)/(cosd(theta(2))*ift(2)) )*1e5;  % [Pa]

Pc_co2br_raw = [1e-5 Pc_co2br_raw];
sg_micp = [0 sg_micp];
Sw_r_clay = 0.3; %0.43;   %0.43
Sg = [0 0.001 0.01:0.01:1-Sw_r_clay];
Pc_co2br  = interp1(sg_micp, Pc_co2br_raw, Sg);
%Pc_co2br_smooth = @(x) 1.162e+7*x.^2 -6.753e+4*x + 3.636e+5; % theta = 0

% Fit power law. Note that we modify the independent term (c) from the best
% fit so that the Pc_e is the same as in the data.
% if theta(1) == 30 
%     Pc_co2br_smooth = @(x) 1.085e+7*x.^2.153 + Pc_co2br(2); % theta = 30
% elseif theta(1) == 70
%     Pc_co2br_smooth = @(x) 4.286e+6*x.^2.153 + Pc_co2br(2); % theta = 70
% elseif theta(1) == 60
%     Pc_co2br_smooth = @(x) 6.265e+6*x.^2.153 + Pc_co2br(2);
% end

pc_fit = fittype('a*exp(b*x)', 'problem', 'a');
[pc_model, gof] = fit(sg_micp', Pc_co2br_raw', pc_fit, 'problem', Pc_co2br_raw(2));
Pc_co2br_smooth = [0; pc_model(Sg(2:end))];

% BC is not a good fit
%pc_e = Pc_co2br(2);
%Sw2 = 1-Sg;
%nSat2 = (Sw2-Sw_r) ./ (1-Sw_r);
%nSat2(end) = 0;
%bc_fit2 = fittype('pc_e*(x)^(-1/lambda)', 'independent', 'x', ...
%                        'coefficients', 'lambda', 'problem', 'pc_e');
%bc_model2 = fit(nSat2(1:end-1)', Pc_co2br(1:end-1)', bc_fit2, 'problem', pc_e);

subplot(1,2,2)
plot(1-Sg, Pc_co2br/1e5, '-', 'color', [0.6 0.6 0.6], 'linewidth', 3, ...
     'DisplayName', 'data');
hold on
plot(1-Sg, Pc_co2br_smooth/1e5, '-r', 'linewidth', 1.5, ...
     'DisplayName', 'power fit');
% plot(1-Sg, [0.01 Pc_co2br_smooth(Sg(2:end))/1e5], '-r', 'linewidth', 1.5, ...
%      'DisplayName', 'power fit');
%plot(1-Sg, bc_model2(nSat2)/1e5, '-r', 'linewidth', 1, ...
%     'DisplayName', 'BC fit');
legend('location', 'southwest')
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$P_\mathrm{c}$ [bar]', latx{:}, 'fontSize', 14)
title('GoM mudstone core $\vert$ z = 3230 m', latx{:})
%set(gca, 'yscale', 'log')
ylim([0 60]), yticks(0:10:60), xlim([0 1])
grid on
xticks(0:.2:1)

%% Relperm 
% 3B.3 test (Gas displacing oil (drainage), Swr > 0)
% Visualize
idB3 = 39:55;
idB3([1 6]) = [];  % too different
sgr = zeros(numel(idB3), 1);
sor = zeros(numel(idB3), 1);
swr = table2array(data(idB3-2, 84));
kr_go_B3 = cell(numel(idB3), 1);
h = figure(2);
subplot(1,2,1)
hold on
for n=1:numel(idB3)
   kr_go_B3{n} = table2array(readtable([pth, fn], 'Sheet', idB3(n))); 
   slen = numel(kr_go_B3{n}(:,1))/2;
   sgr(n) = kr_go_B3{n}(1) + swr(n);
   sor(n) = 1 - kr_go_B3{n}(end, 1) - swr(n);
   kr_go_B3{n}(:, 1) = kr_go_B3{n}(:, 1) + swr(n);
   
   if n==1
   l1 = plot(1-kr_go_B3{n}(slen+1:2*slen, 1), kr_go_B3{n}(slen+1:2*slen, 2), '.-', ...
            'markersize', 4, 'color', [0.8 0.8 0.8], 'DisplayName','3B.3 k$_\mathrm{r,l}$');
   l2 = plot(1-kr_go_B3{n}(1:slen, 1), kr_go_B3{n}(1:slen, 3), '.-', ...
             'markersize', 4, 'color', [0.6 0.6 0.6], 'DisplayName','3B.3 k$_\mathrm{r,g}$');
   %l3 = plot([swr(n) swr(n)], [0, 1], '--', 'color', [0.9 0.9 0.9], ...
   %          'linewidth', 0.5, 'DisplayName', '3B.3 S$_\mathrm{w,i}$');
   else%if n~=6
    plot(1-kr_go_B3{n}(slen+1:2*slen, 1), kr_go_B3{n}(slen+1:2*slen, 2), '.-', ...
         'markersize', 4, 'color', [0.8 0.8 0.8]);
    plot(1-kr_go_B3{n}(1:slen, 1), kr_go_B3{n}(1:slen, 3), '.-', ...
         'markersize', 4, 'color', [0.6 0.6 0.6]);
    %plot([swr(n) swr(n)], [0, 1], '--', 'color', [0.9 0.9 0.9]);
   end
end
latx  = {'Interpreter','latex'};
xlabel('S$_\mathrm{o}$ [-]', latx{:}, 'fontsize', 14)
ylabel('k$_\mathrm{r}$ [-]','fontsize', 14, latx{:})
title('3B.3: g-o (drainage) UCF $\vert$ z = [2452, 2754]m', latx{:})
set(gca,'Xtick',[0 0.2 0.4 0.6 0.8 1])
set(gca,'Ytick',[0 0.25 0.5 0.75 1])
grid on

% Averaged curve
sor_avg = round(mean(sor), 2);
sgr_avg = round(mean(sgr), 2);
sg = sgr_avg:0.005:1-sor_avg;
kr_g_B3_interp = zeros(numel(sg), numel(idB3)+1);
kr_g_B3_interp(:,1) = sg;
kr_o_B3_interp = zeros(numel(sg), numel(idB3)+1);
kr_o_B3_interp(:,1) = sg;
for n=1:numel(kr_go_B3)
    slen = numel(kr_go_B3{n}(:,1))/2;
    kr_g_B3_interp(:,n+1) = interp1(kr_go_B3{n}(1:slen,1), ...
                                    kr_go_B3{n}(1:slen,3), sg);
    kr_o_B3_interp(:,n+1) = interp1(kr_go_B3{n}(slen+1:end,1), ...
                                    kr_go_B3{n}(slen+1:end,2), sg);
end
kr_avg = zeros(numel(sg), 3);
kr_avg(:, 1) = sg;
for j = 1:numel(sg)
    kr_avg(j, 3) = nanmean(kr_g_B3_interp(j, 2:end));
    kr_avg(j, 2) = nanmean(kr_o_B3_interp(j, 2:end));
    if j > 1 && 1-kr_avg(j,1) > 0.53 && kr_avg(j, 2) > kr_avg(j-1, 2)
      kr_avg(j, 2) = 0.99*kr_avg(j-1, 2); 
    end
end
kr_avg = [0.28 1 0.007; 0.29 0.79 0.008; 0.3 0.66 0.009; ...
          0.31 0.58 0.010; 0.32 0.53 0.011; kr_avg];
sgr_avg = 0.28;
sg = [sgr_avg:0.01:round(mean(sgr), 2) ...
      round(mean(sgr), 2)+0.005:0.005:1-sor_avg];
%pc_avg(:, 2) = smooth(pc_avg(:, 1), pc_avg(:, 2), 'rloess');
l3 = plot(1-kr_avg(:,1), kr_avg(:,2), '-', 'color', [0.6 0.6 0.6], ...
          'linewidth', 3, 'DisplayName', 'mean(k$_{ro}$)');
l4 = plot(1-kr_avg(:,1), kr_avg(:,3), '-', 'color', [0.4 0.4 0.4], ...
          'linewidth', 3, 'DisplayName', 'mean(k$_{rg}$)');


% Power fit
So    = 1-sg;
nS    = (So-sor_avg) ./ (1-sor_avg-sgr_avg);
nS(1) = 1; nS(end) = 0;
ofit   = fittype('x^b');  % b~3.4, which gives l~5 if exponent is (2+3l)/l (model of Brooks and Corey, 1966)
o_model = fit(nS', kr_avg(:,2), ofit);
gfit   = fittype('a*x^b');
g_model = fit(1-nS', kr_avg(:,3), gfit);
l5 = plot(So, o_model(nS), '-b', 'linewidth', 1.5, 'DisplayName', 'k$_{og}$ (fit)');
l6 = plot(So, g_model(1-nS), '-r', 'linewidth', 1.5, 'DisplayName', 'k$_{g}$ (fit)');

% GoM synthetic curves: kr, Pc Reservoir
f = fluid;
Sbr   = linspace(f.krPts.og(1, 2), f.krPts.og(1, 3), 50);
Sbs   = linspace(f.krPts.og(2, 2), f.krPts.og(2, 3), 50);
dg = [0.5 0.1 0.1];
l7 = plot(Sbr, f.krOG{1}(Sbr), '-', 'color', 'b', 'linewidth', 0.5, ...
     'DisplayName','k$_\mathrm{r,aq}$');
l8 = plot(Sbr, f.krG{1}(1-Sbr), '-', 'color', 'r', 'linewidth', 0.5, ...
     'DisplayName','k$_\mathrm{r,g}$');
l9 = plot(Sbr, f.krG{4}(1-Sbr), '--', 'color', 'r', 'linewidth', 0.5, ...
     'DisplayName','k$_\mathrm{r,g}^\mathrm{imb}$');
hold off, xlim([0 1]), grid on;
leg = legend([l1, l2, l3, l4, l5, l6, l7, l8, l9] ,latx{:}, 'fontsize', 10);
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
set(h, 'Position', [600, 600, 550, 320])
%set(ax,'color','none')

% kr, Pc Seal
subplot(1,2,2)
plot(Sbs, f.krOG{2}(Sbs), '-', 'color', 'b', 'linewidth', 1);
hold on; plot(Sbs, f.krG{2}(1-Sbs), '-', 'color', 'r', 'linewidth', 1);
hold off,  xlim([0 1]), grid on; ax = gca;
xlabel('S$_\mathrm{aq}$ [-]','fontsize', 15, latx{:})
ylabel('k$_\mathrm{r}$ [-]','fontsize', 15, latx{:})
set(ax,'Xtick',[0 0.2 0.4 0.6 0.8 1]); set(ax,'Ytick',[0 0.25 0.5 0.75 1])
legend('k$_\mathrm{r,aq}$', 'k$_\mathrm{r,g}$', latx{:}, 'fontsize', 10)
title('Caprock')

%% Relperm rescaled to Sgr = 0
% Fit curve re-scaled (same model, but now up to Sw=1)
% For input .DATA file, simply decide what values of Sw you want below, and
% then run code and get krwg and krg. Similar for clay rock with Sw2 below.
%r = (1-min(So))/(max(So)-min(So));
%Sw = min(So) + (So-min(So))*r;
Sw = linspace(min(So),1, 16);
nSw = (Sw-min(Sw))/(1-min(Sw));
krwg = nSw.^o_model.b;
krg = g_model.a*(1-nSw).^g_model.b;

figure(3)
subplot(1,2,1)
hold on
l1 = plot(So, o_model(nS), '-', 'color', [0.8 0.8 0.8], 'linewidth', 1.5, 'DisplayName', 'k$_{og}$ (fit)');
l2 = plot(So, g_model(1-nS), '-', 'color', [0.6 0.6 0.6], 'linewidth', 1.5, 'DisplayName', 'k$_{g}$ (fit)');
l3 = plot(Sw, krwg, '-b', 'linewidth', 1.5, 'DisplayName', 'k$_{wg}$ (scaled)');
l4 = plot(Sw, krg, '-r', 'linewidth', 1.5, 'DisplayName', 'k$_{g}$ (scaled)');
xlabel('S$_\mathrm{w}$ [-]','fontsize', 15, latx{:})
ylabel('k$_\mathrm{r}$ [-]','fontsize', 15, latx{:})
xlim([0 1])
set(gca, 'Xtick',[0 0.2 0.4 0.6 0.8 1]); set(gca,'Ytick',[0 0.25 0.5 0.75 1])
grid on

l5 = plot(Sbr, f.krOG{1}(Sbr), '-', 'color', [0.4 0.4 0.8], 'linewidth', 0.5, ...
     'DisplayName','k$_\mathrm{r,aq}$');
l6 = plot(Sbr, f.krG{1}(1-Sbr), '-', 'color', [0.8 0.4 0.4], 'linewidth', 0.5, ...
     'DisplayName','k$_\mathrm{r,g}$');
l7 = plot(Sbr, f.krG{4}(1-Sbr), '--', 'color', [0.8 0.4 0.4], 'linewidth', 0.5, ...
     'DisplayName','k$_\mathrm{r,g}^\mathrm{imb}$');
leg = legend([l1, l2, l3, l4, l5, l6, l7] ,latx{:}, 'fontsize', 10, 'location', 'northwest');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
title('Sand-rich material')

% Clay curves with same exponent for water
%Sw2 = linspace(f.krPts.og(2, 2),1, 16);
Sw2 = linspace(Sw_r_clay,1, 16);
nSw2 = (Sw2-min(Sw2))/(1-min(Sw2));
krwg2 = nSw2.^o_model.b;
krg2 = (g_model.a-0.1)*(1-nSw2).^g_model.b;

subplot(1,2,2)
hold on
%l1 = plot(Sw2, krwg2, '-b', 'linewidth', 1.5, 'DisplayName', 'k$_{wg}$ (scaled)');
%l2 = plot(Sw2, krg2, '-r', 'linewidth', 1.5, 'DisplayName', 'k$_{g}$ (scaled)');
l3 = plot(Sbs, f.krOG{2}(Sbs), '-', 'color', [0.4 0.4 0.8], 'linewidth', 0.5, 'DisplayName','k$_\mathrm{r,aq}$');
l4 = plot(Sbs, f.krG{2}(1-Sbs), '-', 'color', [0.8 0.4 0.4], 'linewidth', 0.5, 'DisplayName','k$_\mathrm{r,g}$');
hold off,  xlim([0 1]), grid on; ax = gca;
xlabel('S$_\mathrm{aq}$ [-]','fontsize', 15, latx{:})
ylabel('k$_\mathrm{r}$ [-]','fontsize', 15, latx{:})
set(ax,'Xtick',[0 0.2 0.4 0.6 0.8 1]); set(ax,'Ytick',[0 0.25 0.5 0.75 1])
title('Clay-rich material')
leg = legend([l3, l4] ,latx{:}, 'fontsize', 10, 'location', 'northwest');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
set(h, 'Position', [600, 600, 550, 320])

%% Comparison with CO2-brine measurements in literature
% We use capillary-dominated experiments to compare, because XOM data were
% obtained for VL conditions (high flow rates). However, most parts of the
% simulation model will be dominated by capillarity/gravity.

% Reynolds & Krevor, WRR (2015)
% Bentheimer sst, with about 95% Q, 4% F, 1% clay (scCO2/brine)
rey15CL = {[0.7782 0.7185 0.6582 0.6084 0.5419 0.4544 0.3525;    % Sw, exp2
            0.0015 0.0043 0.0137 0.0456 0.0864 0.1138 0.1806;    % krg
            0.3139 0.1885 0.1342 0.1002 0.0563 0.0155 0.0048]    % krw
           [0.8293 0.8107 0.8037 0.7589 0.7217 0.6746 0.6506 0.5569 0.4673 0.3649; % exp 4
            0.0000 0.0001 0.0005 0.0012 0.0034 0.0087 0.0237 0.0575 0.0839 0.1158;
            0.4755 0.4604 0.3630 0.2500 0.2008 0.1550 0.1277 0.0690 0.0216 0.0022]
           [0.8969 0.8762 0.8434 0.7804 0.7426 0.7055 0.6158 0.5576 0.4873 0.3349 % exp 6
            0.0006 0.0015 0.0030 0.0055 0.0101 0.0195 0.0603 0.0905 0.1311 0.2200
            0.7056 0.5805 0.4577 0.3274 0.2435 0.1914 0.0980 0.0589 0.0246 0.0007]};

% Akbarabadi & Piri, AWR (2013)
% Berea sst, about 10% clay (scCO2/brine, steady-state)
akb13tr = {[1 0.997 9.995 0.988 0.986 0.969 0.928 0.836 0.77 0.724 0.676 ...
            0.644 0.621 0.577 0.553 0.515;    % Sw, exp 24
            1 0.484 0.420 0.404 0.312 0.266 0.155 0.114 0.093 0.051 0.037 ...
            0.022 0.020 0.007 0.007 0;        % krw
            0 0.005 0.006 0.011 0.012 0.021 0.025 0.037 0.056 0.062 0.09 ...
            0.104 0.128 0.141 0.167 0.187]    % krg
           [1 0.995 0.994 0.992 0.988 0.963 0.951 0.853 0.820 0.760 0.731 ...
            0.692 0.664 0.633 0.595 0.575 0.535;    % exp 25
            1 0.494 0.484 0.461 0.347 0.246 0.167 0.129 0.115 0.086 0.052 ...
            0.041 0.023 0.018 0.008 0.007 0;
            0 0.005 0.007 0.013 0.013 0.019 0.027 0.041 0.046 0.063 0.076 ...
            0.1 0.112 0.119 0.153 0.166 0.191]};

% Pini & Benson, WRR (2013), scCO2/water
% Berea sst , Fig 8 (BC fit not good, use datapoints)
% pb13.swirr = 0.325;
% pb13.lambda = 2.36;
% pb13.Sw = pb13.swirr:.01:1;
% pb13.nSw = (pb13.Sw - pb13.swirr) / (1 - pb13.swirr);
% pb13.krn = (1-pb13.nSw).^2.*(1-pb13.nSw.^((2+pb13.lambda)/pb13.lambda));
pb13CL = [0.64 0.63 0.6 0.54 0.5 0.46 0.43 0.4 0.39 0.38;
          0.07 0.08 0.11 0.18 0.26 0.32 0.42 0.58 0.72 0.91];

h = figure(4);
hold on
l1 = plot(rey15CL{1}(1,:), rey15CL{1}(2,:), '-s', ...
          'color', [0.8 0.4 0.4], 'DisplayName', 'RK15, Bentheimer');  
plot(rey15CL{2}(1,:), rey15CL{2}(2,:), '-s', 'color', [0.8 0.4 0.4]);
plot(rey15CL{3}(1,:), rey15CL{3}(2,:), '-s', 'color', [0.8 0.4 0.4]);
plot(rey15CL{1}(1,:), rey15CL{1}(3,:), '-s', 'color', [0.4 0.4 0.8])
plot(rey15CL{2}(1,:), rey15CL{2}(3,:), '-s', 'color', [0.4 0.4 0.8])
plot(rey15CL{3}(1,:), rey15CL{3}(3,:), '-s', 'color', [0.4 0.4 0.8])
l2 = plot(akb13tr{1}(1,:), akb13tr{1}(3,:), '--or', 'markerSize', 4, ...
          'DisplayName', 'AP13, Berea', 'color', [0.8 0.4 0.4]);
plot(akb13tr{2}(1,:), akb13tr{2}(3,:), '--or', 'markerSize', 4, 'color', [0.8 0.4 0.4])
plot(akb13tr{1}(1,:), akb13tr{1}(2,:), '--ob', 'markerSize', 4, 'color', [0.4 0.4 0.8])
plot(akb13tr{2}(1,:), akb13tr{2}(2,:), '--ob', 'markerSize', 4, 'color', [0.4 0.4 0.8])
l3 = plot(pb13CL(1,:), pb13CL(2,:), '-.d', 'markerSize', 6, ...
          'DisplayName', 'PB13, Berea', 'color', [0.8 0.4 0.4]);
l4 = plot(Sw, krwg, '-b', 'linewidth', 1.5, 'DisplayName', 'k$_{wg}$ (scaled)');
l5 = plot(Sw, krg, '-r', 'linewidth', 1.5, 'DisplayName', 'k$_{g}$ (scaled)');
hold off
xlabel('S$_\mathrm{w}$ [-]','fontsize', 15, latx{:})
ylabel('k$_\mathrm{r}$ [-]','fontsize', 15, latx{:})
xlim([0 1])
set(gca, 'Xtick',[0 0.2 0.4 0.6 0.8 1]); set(gca,'Ytick',[0 0.25 0.5 0.75 1])
grid on
leg = legend([l1, l2, l3, l4, l5] ,latx{:}, 'fontsize', 10, 'location', 'northwest');
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
set(h, 'Position', [200, 200, 400, 400])


%% Re-scale models to krg max = 1
% add extra points at new residual saturation first, then fit.
% Sand
Sw_r_sand = 0.16;
Sw_new = [Sw_r_sand 0.2 Sw];
Sg_new = 1 - Sw_new;
krg_new = [1.05 0.8 krg];
nSw = (Sw_new - min(Sw_new)) ./ (1 - min(Sw_new));
gfit   = fittype('a*x^b');
g_model = fit(1-nSw', krg_new', gfit);

Sw_new2 = [Sw_new(1) Sw(2:end)];
krw_new = [0 krwg(2:end)];
nSw2 = (Sw_new2 - min(Sw_new2)) ./ (1 - min(Sw_new2));
wfit   = fittype('x^b');
w_model = fit(nSw2', krw_new', wfit);

sw_pc = [pc_avg(1:end-1,1); 1.001*Sw_r_sand];
bcv = [bc_vals; 1e6];
nSw_pc = (sw_pc-Sw_r_sand) ./ (1-Sw_r_sand);
pc_t = bc_vals(1);
bc_fit = fittype('pc_t*(x)^(-1/lambda)', 'independent', 'x', ...
                 'coefficients', 'lambda', 'problem', 'pc_t');
bc_model2 = fit(nSw_pc, bcv, bc_fit, 'problem', pc_t);

% Clay
Swc_new = linspace(Sw_r_clay,1, 32);
nSwc_new = (Swc_new-min(Swc_new))/(1-min(Swc_new));

% Figs
figure(3); 
subplot(1,2,1); 
hold on; plot(Sw_new2, w_model(nSw2), '-c', 'displayName', '$k_\mathrm{r,aq}$ (fin)', ...
              'linewidth', 1.25);
plot(Sw_new, g_model(1-nSw), '-', 'color', [247, 129, 10]/255, ...
     'displayName', '$k_\mathrm{r,g}$ (fin)', 'linewidth', 1.25);
hold off
xlim([0 1])
subplot(1,2,2)
hold on; plot(Swc_new, w_model(nSwc_new), '-c', 'displayName', '$k_\mathrm{r,aq}$ (fin)', ...
              'linewidth', 1.25);
plot(Swc_new, g_model(1-nSwc_new), '-', 'color', [247, 129, 10]/255, ...
     'displayName', '$k_\mathrm{r,g}$ (fin)', 'linewidth', 1.25);

figure(18); hold on; 
plot(sw_pc, bcv, 'sk')
swp = [0.16005:0.00005:0.1601 0.162 0.165 0.17 0.18:0.02:1];
nswp = (swp-0.16) ./ (1-0.16);
plot(swp, bc_model2(nswp)); 
set(gca,'YScale', 'log')

% caprock kr: same model as sand, rescaled to residual gas sat. Increase
% sgmax to 0.3 or so. Use data for Pc curve.


%% Write data to file
% sgmax = [1-swr_avg_sand, 1-Sw_r_clay];
sgmax = [1-Sw_r_sand, 1-Sw_r_clay];
sg_sand = linspace(0, sgmax(1), 30);
sg_sand = [sg_sand(1) 0.001 sg_sand(2:end-1) sgmax(1)-0.001 sg_sand(end)];
sg_clay = linspace(0, sgmax(2), 30);
sg_clay = [sg_clay(1) 0.001 sg_clay(2:end-1) sgmax(2)-0.001 sg_clay(end)];

% Sand Pc
%nSat_sand = ((1-sg_sand)-swr_avg_sand) ./ (1-swr_avg_sand-sg_sand(1));
%pc_sand = [0; bc_model(nSat_sand(1:end-1))/1e5];
nSat_sand = ((1-sg_sand)-Sw_r_sand) ./ (1-Sw_r_sand);
pc_sand = [0; bc_model2(nSat_sand(2:end-1))/1e5; 10];

% Clay Pc
%pc_clay = [0 Pc_co2br_smooth(sg_clay(2:end))/1e5];
pc_clay = [0; pc_model(sg_clay(2:end))/1e5];

% Sand kr
% Sw1    = 1-sg_sand;
% nSwg = (Sw1-min(Sw))/(1-min(Sw));   % here we want to extend the kr curve to the sgmax from Pc sand (not rescale) 
% nSw1 = (Sw1-min(Sw1))/(1-min(Sw1));
% kr_w = nSw1.^o_model.b;
% kr_n = g_model.a*(1-nSwg).^g_model.b;
kr_w = w_model(nSat_sand);
kr_w(end) = 0;
kr_n = g_model(1-nSat_sand);

% Clay kr
% Sw2 = 1-sg_clay;
% nSw2 = (Sw2-min(Sw2))/(1-min(Sw2));
% kr_w2 = nSw2.^o_model.b;
% kr_n2 = (g_model.a-0.1)*(1-nSw2).^g_model.b;
nSat_clay = ((1-sg_clay)-Sw_r_clay) ./ (1-Sw_r_clay);
kr_w2 = w_model(nSat_clay); 
kr_w2(end) = 0;
kr_n2 = g_model(1-nSat_clay);

% plots
figure(5)
subplot(1,2,1)
hold on
plot(1-sg_sand, kr_w, '-b', 'linewidth', 1.25)
plot(1-sg_sand, kr_n, '-r', 'linewidth', 1.25)
hold off
grid on
xlim([0 1]); ylim([0 1])
xticks(0:.2:1); yticks(0:.2:1);
xlabel('$S_\mathrm{w}$ [-]', 'interpreter', 'latex', 'fontSize', 12)
ylabel('$k_\mathrm{r}$ [-]', 'interpreter', 'latex', 'fontSize', 12)
legend('brine', 'CO_2-rich')
subplot(1,2,2)
plot(1-sg_sand, pc_sand, '-ok', 'linewidth', 1.25, 'markersize', 3)
grid on
xlim([0 1]); ylim([10^-2 10^1]); set(gca, 'YScale', 'log')
xticks(0:.2:1);
xlabel('$S_\mathrm{w}$ [-]', 'interpreter', 'latex', 'fontSize', 12)
ylabel('$P_\mathrm{c}$ [bar]', 'interpreter', 'latex', 'fontSize', 12)

figure(6)
subplot(1,2,1)
hold on
plot(1-sg_clay, kr_w2, '-b', 'linewidth', 1.25)
plot(1-sg_clay, kr_n2, '-r', 'linewidth', 1.25)
hold off
grid on
xlim([0 1]); ylim([0 1])
xticks(0:.2:1); yticks(0:.2:1);
xlabel('$S_\mathrm{w}$ [-]', 'interpreter', 'latex', 'fontSize', 12)
ylabel('$k_\mathrm{r}$ [-]', 'interpreter', 'latex', 'fontSize', 12)
%legend('brine', 'CO_2-rich')
subplot(1,2,2)
plot(1-sg_clay, pc_clay, '-ok', 'linewidth', 1.25, 'markersize', 3)
grid on
xlim([0 1]); %ylim([10^-2 10^1]); set(gca, 'YScale', 'log')
xticks(0:.2:1);
xlabel('$S_\mathrm{w}$ [-]', 'interpreter', 'latex', 'fontSize', 12)
ylabel('$P_\mathrm{c}$ [bar]', 'interpreter', 'latex', 'fontSize', 12)

figure(7)
vpar_fit = [4.4 2.8 0.18];  % from upscaling
sg_ups = 0:.01:(1-vpar_fit(3));
sw_ups = 1-sg_ups;
sw_ups_n = (sw_ups - min(sw_ups))/(1 - min(sw_ups));
sg_ups_n = 1 - sw_ups_n;
krw_ups = sw_ups_n.^vpar_fit(1);
krg_ups = sg_ups_n.^vpar_fit(2);
hold on
p1 = plot(1-sg_sand, kr_w, '-b', 'linewidth', 1, 'DisplayName', '$k_\mathrm{rw}$ (sand)');
p2 = plot(1-sg_sand, kr_n, '-r', 'linewidth', 1, 'DisplayName', '$k_\mathrm{rg}$ (sand)');
p3 = plot(1-sg_clay, kr_w2, '--b', 'linewidth', 1, 'DisplayName', '$k_\mathrm{rw}$ (clay)');
plot(1-sg_clay, kr_n2, '--r', 'linewidth', 1)
p4 = plot(sw_ups, krw_ups, '.b', 'linewidth', 1, 'DisplayName', '$k_\mathrm{rw}$ (ups)');
plot(sw_ups, krg_ups, '.r', 'linewidth', 1)
hold off
grid on
xlim([0 1]); ylim([0 1])
xticks(0:.2:1); yticks(0:.2:1);
xlabel('$S_\mathrm{w}$ [-]', 'interpreter', 'latex', 'fontSize', 14)
ylabel('$k_\mathrm{r}$ [-]', 'interpreter', 'latex', 'fontSize', 14)
legend([p1 p2 p3 p4],'interpreter','latex','fontsize', 11)

% table
fileName  = 'SGOF_theta30_rescaledTo1_noPc';
satg      = [sg_sand NaN sg_clay]';
krn       = [kr_n; NaN; kr_n2];
krw       = [kr_w; NaN; kr_w2];
%Pc        = [pc_sand; NaN; pc_clay];
Pc        = [zeros(numel(pc_sand), 1); NaN; zeros(numel(pc_clay), 1)];
        
tab = table(satg, krn, krw, Pc);
tab.Properties.VariableNames = {'SGAS' 'KRG' 'KRWG' 'PCOG'};
writetable(tab, [fileName '.txt'], 'Delimiter', 'tab');
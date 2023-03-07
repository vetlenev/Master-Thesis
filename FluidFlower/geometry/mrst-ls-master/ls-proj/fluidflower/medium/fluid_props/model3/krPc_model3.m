%% Model 1 relperm and Pc
%
% Summary:
%
% Description:
%   Pc: 
%   krw (drainage):
%   krg (drainage):
%   krg (imbibition): 
%s
%      
%
clear, close all
model_case = 3;
sg_pce = 1e-4;  
sandtype = 'C';
pcmult = [2 1.7 1 1];     % [ESF C E F]
np = 16;                % number of points in output

%% Sand data (see getRock.m) 
d_mm = [0.1, 0.66, 1.45, 1.77];     % ESF, C, E, F
porov = [0.435, 0.435, .45, .44];   % "
permv = [44, 473, 2005, 4259];      % "

switch sandtype
    case 'ESF'
        ids = 1;
    case 'C'
        ids = 2;
    case 'E'
        ids = 3;
    case 'F'
        ids = 4;
end

%% Init properties and residual water saturation
% We fit Pc and kr using the data (including swc from data). After fitting,
% we re-scale the Pc and kr curves to the new saturation space with Swc for 
% our sands.
k = permv(ids);
n = porov(ids);
% modif krel gas and krel wat from data in description.pdf
swcv = [0.32 0.14 0.12 0.12];
sgtv = [0.14 0.10 0.06 0.13];
pcev = [15, 3, 0.33, 0].*pcmult;   %   mbar (E sand based on trend (not reported))
swc = swcv(ids);
sgt = sgtv(ids);
pce = pcev(ids);


%% Pc
% BC fit to approximate avg of data in Plug & Bruining, AWR (2007; Fig. 4).
sandpack = 'coarse';
if strcmp(sandpack, 'fine')
    % Avg T = 22 C, Fine sand (160 < D_50 < 210 micrometers), 
    % avg sandpack porosity 0.35. Fig. 4
    data.pc_vals_bar = fliplr([28 30 40 50 60 70 80 90 100])/1000;
    data.pc_sw = fliplr([1 0.98 0.71 0.48 0.38 0.34 0.31 0.28 0.26]);
elseif strcmp(sandpack, 'coarse')
    % Avg T = 25 C, Coarse sand (360 < D_50 < 410 micrometers), 
    % avg sandpack porosity 0.37, perm 2e-10 m^2. Fig. 5.
    data.pc_d50 = mean([0.360 0.410]);
    data.pc_n = 0.37;
    data.pc_k = 2e-10/darcy;
    data.pc_vals_bar = fliplr([11 12.5 14 16 17 18 18.5 19 20.5 23 32 45])/1000;
    data.pc_sw = fliplr([1 0.97 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.07]);
end
s_sub = 2e-2;        % to use the first point we need it to be > swc 
swc_pcdat = min(data.pc_sw)-s_sub;

% Fit Brooks-Corey (1964) model
% https://petrowiki.spe.org/Capillary_pressure_models
pc_t          = min(data.pc_vals_bar);                      % threshold pc
nSat          = (data.pc_sw-swc_pcdat) ./ (1-swc_pcdat);
bc_fit        = fittype('pc_t*(x)^(-1/lambda)', 'independent', 'x', ...
                        'coefficients', 'lambda', 'problem', 'pc_t');
[bc_model, gof_pc] = fit(nSat(2:end)', data.pc_vals_bar(2:end)', bc_fit, 'problem', pc_t);

sw = linspace(swc_pcdat+s_sub, 1, 50);
ns = (sw-swc_pcdat) ./ (1-swc_pcdat);
pc_fit = bc_model(ns);

% Figure
h = figure(1);
subplot(1,2,1)
l1 = plot(data.pc_sw, data.pc_vals_bar*1000, 'sk', ...
         'DisplayName', 'Avg Plug \& Bruining (2007)');
hold on
l2 = plot(sw, pc_fit*1000, '-k', 'linewidth', 1.5, ...
         'DisplayName', 'BC64 fit');
latx  = {'Interpreter','latex'};
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$P_\mathrm{c}$ [mbar]', latx{:}, 'fontSize', 14)
%set(gca, 'yscale', 'log')
ylim([0 60]), yticks(0:5:100)
xlim([0 1]), xticks(0:.1:1)
grid on
leg = legend([l1, l2], 'location', 'northeast', latx{:});
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.8]));
hold off
title(['$\overline{d_{50}}$ = ' num2str(data.pc_d50), 'mm $\mid$ ' ...
       '$\phi$ = ' num2str(data.pc_n), ...
       ' $\mid$ $k$ = ', num2str(data.pc_k), 'D'], latx{:});


%% Relperm Drainage
% DiCarlo et al., SPE (2000) Fig.4 and 5
data.krg = [.025 .04 .07 .04 .075 .06 .125 .24 .16 .18 .255 .35 .44 ...
            .52 .62 .64 .81 .8 .875 .87 .95 .82 .99];          
data.krg_sg = [.11 .12 .12 .14 .14 .17 .245 .28 .31 .34 .34 .375 .435 ...
            .5 .6 .62 .72 .75 .79 .81 .825 .85 .855];
swc_krdat = 1-0.86;
nsw_krdat = ((1-data.krg_sg)-swc_krdat)/(1-swc_krdat);
g_fit = fittype('a*x^b', 'independent', 'x');
g_model = fit(1-nsw_krdat', data.krg', g_fit);


sw = swc_krdat:.01:1;
nsw = (sw-swc_krdat)/(1-swc_krdat);
krg = g_model(1-nsw);
krw = sw.^5; 
w_fit = fittype('x^b', 'independent', 'x');
w_model = fit(nsw', krw', w_fit);


subplot(1,2,2)
l1 = plot(1-data.krg_sg, data.krg, 'sk', ...
         'DisplayName', 'DiCarlo et al. (2000)');
hold on
l2 = plot(sw, krg, '-r', 'linewidth', 1.5, ...
         'DisplayName', '$0.97(1-S_w^*)^{1.39}$');
l3 = plot(sw, krw, '-b', 'linewidth', 1.5, ...
         'DisplayName', '$S_w^5$ (DiCarlo et al.)');
l4 = plot(sw, w_model(nsw), '--c', 'linewidth', 1, ...
         'DisplayName', '$(S_w^*)^{4.192}$');
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$k_r$ [-]', latx{:}, 'fontSize', 14)
%set(gca, 'yscale', 'log')
ylim([0 1]), yticks(0:.1:1)
xlim([0 1]), xticks(0:.1:1)
grid on
leg = legend([l1, l2, l3, l4], 'location', 'northeast', latx{:});
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.8]));
hold off
set(h, 'Position', [200, 200, 650, 325])


%% Rescale Pc and Kr

% Swc
sw = zeros(np, 1);
sw([1:end-2 end]) = linspace(swc, 1, np-1);
sw(end-1) = 1-sg_pce;         % Has to be > 0.999 otherwise very low sat may enter the seal.
ns = (sw-swc) ./ (1-swc);
nspcmin = 0.01;
nspcfit = [nspcmin; ns(2:end)];
pc_fit_scaled = [bc_model(nspcmin); bc_model(ns(2:end))];
krg_fit_scaled = g_model(1-ns);
krw_fit_scaled = w_model(ns);
pc_fit_scaled = pc_fit_scaled*pce/(pc_fit_scaled(end-1)*1000);


% Figure
h = figure(2);
subplot(1,2,1)
l1 = plot(sw, pc_fit_scaled*1000, 'ok', 'markersize', 6, ...
         'DisplayName', '$P_{c}$: BC64 fit, Swc');
latx  = {'Interpreter','latex'};
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$P_\mathrm{c}$ [mbar]', latx{:}, 'fontSize', 14)
%set(gca, 'yscale', 'log')
ylimv = round(max(pc_fit_scaled), 3)*1000 + 5;
switch sandtype
    case 'ESF', spac = 20;
    case 'C', spac=5;
    case 'D', spac=5;
    case 'E', spac=5;
    case 'F', spac=3;
    case 'G', spac=2;
end
ylim([0 ylimv]), yticks(0:spac:ylimv)
xlim([0 1]), xticks(0:.1:1)
grid on
leg = legend(l1, 'location', 'northeast', latx{:});
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.8]));
hold off

subplot(1,2,2)
krg_name = '$k_{r\mathrm{g}}$: $0.97(1-S_w^*)^{1.39}$, Swc';
krw_name = '$k_{r\mathrm{w}}$: $(S_w^*)^{4.192}$, Swc';
l1 = plot(sw, krg_fit_scaled, '-r', 'linewidth', 2, ...
         'DisplayName', krg_name);
hold on
l2 = plot(sw, krw_fit_scaled, '-b', 'linewidth', 2, ...
         'DisplayName', krw_name);
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$k_r$ [-]', latx{:}, 'fontSize', 14)
%set(gca, 'yscale', 'log')
ylim([0 1]), yticks(0:.1:1)
xlim([0 1]), xticks(0:.1:1)
grid on
leg = legend([l1, l2], 'location', 'northeast', latx{:});
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.8]));
hold off
set(h, 'Position', [200, 200, 650, 325])


%% Relperm Imbibition
% Calculate trapped gas saturation
% Data from Pentland et al., SPE (2010) for experiments in sandpacks. Land
% model is not the one that fits best but it makes sense to start here,
% also because it is consistent with Pc curve from Brooks Corey with
% parameter lambda that is also used to calculate krib.
data.sgi = [3 4 6 6 6.5 7 7 14 15.5 26 31.125 38.7 40 42 ...
            44 46.25 48 50 52 58 58.5 59 59 59.5 60.5 63 64 65 ...
            66 68 70 71 72 72.3 73.75 73.4 76.25 77.5 79.5]./100;
data.sgt = [0.7 1 1.3 2 1.6 1.5 3.2 3 6.1 7.5 9.3 11 9.6 10.5...
            11.7 13.1 12.2 12 12.3 11.6 12.4 13.1 13.3 12.8 12.5 13.7 12.4 13.3 ...
            11.9 11.4 13.2 12.6 13.1 12.7 13.1 13.9 13 12.7 12.8]./100;
data_swc = 1 - round(max(data.sgi),1);
data_nsgi = data.sgi ./ (1-data_swc);      % same as 1 - ((sw-swc)/(1-swc))
data_nsgt = data.sgt ./ (1-data_swc);
% Here, model result is also normalized (saturation) so it needs to be 
% de-normalized below!
sgt_fit   = fittype('x/(1+a*x)', 'independent', 'x');
[sgt_model, gof] = fit(data_nsgi', data_nsgt', sgt_fit); % C=5.2, R^2 = 0.89
sgi = 1 - sw;
sgt_vals = sgt_model(1-ns)*(1-data_swc);  % de-normalize result.

h = figure(3);
hold on
l1 = plot(data.sgi, data.sgt, 'sk', 'DisplayName', 'Pentland et al. (2010)');
l2 = plot(sgi, sgt_vals, '--k', 'DisplayName', ...
          '$S_\mathrm{gt}^* = S_\mathrm{gi}^*/(1+5.2S_\mathrm{gi}^*)$');
grid on
hold off
xlabel('S_{gi}')
ylabel('S_{gt}')
xlim([0 1]), xticks(0:.2:1)
ylim([0 0.15]), yticks(0:0.025:0.15)
legend([l1, l2], 'location', 'southeast', latx{:})
set(h, 'Position', [200, 200, 350, 275])

% Calculate imbibition curves
% Gas (Land, 1968)
sgi_max = max(sgi);
nsgi_max = sgi_max / (1-swc); % this is always 1 for bounding ib curve
if model_case == 1
    nsgt_max = sgt_model(nsgi_max);
    sgt_max = nsgt_max*(1-swc);
elseif model_case == 3
    sgt_max = sgt;
    nsgt_max = sgt_max/(1-swc);
end
%nsgt_max = 0.4;                              % Land (1968) Fig. 2
C = 1/nsgt_max - 1;                           % Land (1968) Eq. 2
sgib = linspace(sgt_max, sgi_max, np);
nsgib = (sgib)/(1-swc);            
nsgibF = 0.5*((nsgib-nsgib(1)) + sqrt((nsgib-nsgib(1)).^2 + ...
                                      (4/C)*(nsgib-nsgib(1))));  % Land (1968) Eq. 4
e = 2/bc_model.lambda + 3;
%e = 4;                                      % Land (1968) Fig. 2
krgib = nsgibF.^2.*(1-(1-nsgibF).^(e-2));
krgib = max(krg_fit_scaled)*krgib;             % our krg is not 1 at sgi;

% Water and Pcog (no hysteresis)
nswib = 1 - nsgib;
krwib = w_model(nswib);
nspcfitib = [nswib(1:end-1) nspcmin];
pc_fit_scaled_ib = interp1(nspcfit, pc_fit_scaled, nspcfitib)';

figure(2)
subplot(1,2,1)
hold on
plot(1-sgib, pc_fit_scaled_ib*1e3, 'sk', 'markerfacecolor', 'k', ...
     'markerSize', 4, 'DisplayName', '$P_{c}^\mathrm{ib}$');
hold off
subplot(1,2,2)
hold on
plot(1-sgib, krgib, '--r', 'linewidth', 1, 'DisplayName', '$k_{r\mathrm{g}}^\mathrm{ib}$: Land (1968), Eq. 5')
plot(1-sgib, krwib, '--c', 'linewidth', 1, 'DisplayName', '$k_{r\mathrm{w}}^\mathrm{ib}$')
sgtitle(['Sand ' sandtype ', $\overline{d}$ = ' ...
         num2str(d_mm(ids)) 'mm'], latx{:})
hold off


%% Generate tables for .DATA files
hdrs = {'SGAS', 'KRG', 'KROG', 'PCOG'};
vartyp  = cellstr(repmat('double', numel(hdrs), 1));
t.drainage = table('Size', [np, 4], 'VariableTypes', vartyp, ...
                   'VariableNames', hdrs); 
t.imbibition = table('Size', [np, 4], 'VariableTypes', vartyp, ...
                   'VariableNames', hdrs);
t.drainage = flipud([(1-sw), krg_fit_scaled, krw_fit_scaled, ...
                     [pc_fit_scaled(1:end-1); 0]]);
t.imbibition = [sgib', krgib', krwib, pc_fit_scaled_ib];
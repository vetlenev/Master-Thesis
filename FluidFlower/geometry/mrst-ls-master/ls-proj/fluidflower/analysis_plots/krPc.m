%% Model 1 relperm and Pc
%
% Plot analysis
%
clear, close all
model_case = 2;        % 1 medium match, 2 model 2, 32 model 3
sg_pce = 1e-4;  
pcmult = 1;
np = 2^5;               % number of points in output
plot_data = false;

%% Sand data (see getRock.m)
d_mm = [0.2, 0.66, 1.05, 1.45, 1.77, 2.51];
%d_mm = [0.1, 0.66, 1.05, 1.45, 1.77, 3.0];  % ESF, C, D, E, F, G [mm]
if model_case == 1
    porov = [.37, .38, .39, .39, .39, .42];
    perm_model = @(x,y) 1.225e+4*x.^2.*y.^3;
    permv = perm_model(d_mm, porov);            % D
else % data in description.pdf
    porov = [0.435, 0.435, .44, .45, .44, .45];
    permv = [44, 473, 1110, 2005, 4259, 9580];  % D
end

%% Init properties and residual water saturation
% We fit Pc and kr using the data (including swc from data). After fitting,
% we re-scale the Pc and kr curves to the new saturation space with Swc for 
% our sands.
k = permv;
phi = porov;
if model_case == 1 || model_case == 2
    % Swc is computed from Timur, The Leadig Edge (1968)
    swc = round(0.01*(3.5*((phi*1e2).^1.26)./((k*1e3).^0.35)-1), 2); % in [0 1]
elseif model_case == 3 % data in description.pdf
    swcv = [0.32 0.14 0.12 0.12 0.12 0.10];
    krg_swcv = [0.09 0.05 0.02 0.1 0.11 0.16];
    sgtv = [0.14 0.10 0.08 0.06 0.13 0.06];
    krw_sgtv = [0.71 0.93 0.95 0.93 0.72 0.75];
    pcev = [15, 3, 1, 0.33, 0, 0];   %   mbar (E sand based on trend (not reported))
    swc = swcv(ids);
    sgt = sgtv(ids);
    krg_swc = krg_swcv(ids);
    krw_sgt = krw_sgtv(ids);
    pce = pcev(ids);
elseif model_case == 32 % modif krel gas and krel wat from data in description.pdf
    swcv = [0.32 0.14 0.12 0.12 0.12 0.10];     % model 3
    sgtv = [0.14 0.10 0.08 0.06 0.13 0.06];
%     swcv = [0.32 0.14 0.02 0.12 0.12 0.01];   % model 2 (only D, G)
%     sgtv = [0.14 0.10 0.1573 0.06 0.13 0.16];
    pcev = [15, 3, 1, 0.33, 0, 0]*pcmult;   %   mbar (E sand based on trend (not reported))
    swc = swcv;
    sgt = sgtv;
    pce = pcev;
end


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
if plot_data
    h = figure(1);
    subplot(1,5,1)
    l1 = plot(data.pc_sw, data.pc_vals_bar*1000, 'sk', ...
             'DisplayName', 'Plug \& Bruining (2007)');
    hold on
    l2 = plot(sw, pc_fit*1000, '-k', 'linewidth', 1.5, ...
             'DisplayName', 'BC64 fit');
    latx  = {'Interpreter','latex'};
    xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
    ylabel('$p_\mathrm{c}$ [mbar]', latx{:}, 'fontSize', 14)
    %set(gca, 'yscale', 'log')
    ylim([0 60]), yticks(0:5:100)
    xlim([0 1]), xticks(0:.2:1)
    grid on
    leg = legend([l1, l2], 'location', 'northeast', latx{:});
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.8]));
    hold off
    % title(['$\overline{d_{50}}$ = ' num2str(data.pc_d50), 'mm $\mid$ ' ...
    %        '$\phi$ = ' num2str(data.pc_n), ...
    %        ' $\mid$ $k$ = ', num2str(data.pc_k), 'D'], latx{:});
    %title('p_c fit', 'fontweight', 'normal');
end


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


if plot_data
    subplot(1,5,3)
    l1 = plot(1-data.krg_sg, data.krg, 'sk', ...
             'DisplayName', 'DiCarlo et al. (2000)');
    hold on
    l2 = plot(sw, krg, '-r', 'linewidth', 1.5, ...
             'DisplayName', '$0.97(1-S_w^*)^{1.39}$');
    l3 = plot(sw, krw, '-b', 'linewidth', 1.5, ...
             'DisplayName', '$S_w^5$');
    l4 = plot(sw, w_model(nsw), '--c', 'linewidth', 1, ...
             'DisplayName', '$(S_w^*)^{4.192}$');
    xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
    ylabel('$k_r$ [-]', latx{:}, 'fontSize', 14)
    %set(gca, 'yscale', 'log')
    ylim([0 1]), yticks(0:.1:1)
    xlim([0 1]), xticks(0:.2:1)
    grid on
    leg = legend([l1, l2, l3, l4], 'location', 'northeast', latx{:});
    set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.8]));
    hold off
    %title('k_r fit', 'fontweight', 'normal');
end


%% Rescale Pc and Kr
m = numel(k);
sw = zeros(np, m);
for n=1:m
% Swc
sw([1:end-2 end], n) = linspace(swc(n), 1, np-1);
sw(end-1, n) = 1-sg_pce;         % Has to be > 0.999 otherwise very low sat may enter the seal.
ns = (sw(:,n)-swc(n)) ./ (1-swc(n));
nspcmin = 0.01;
nspcfit = [nspcmin; ns(2:end)];
pc_fit_scaled1 = [bc_model(nspcmin); bc_model(ns(2:end))];
krg_fit_scaled(:,n) = g_model(1-ns);
krw_fit_scaled(:,n) = w_model(ns);

if model_case == 1 || model_case == 2
    % Use Leverett-J scaling to re-scale curve to new poro-perm
    pc_fit_scaled(:,n) = pcmult*pc_fit_scaled1.*sqrt((data.pc_k*phi(n))./(data.pc_n*k(n)));
elseif model_case == 3
    % Scale Pc and krg based on Pce and endpoint (swc) krg, respectively
    pc_fit_scaled = pc_fit_scaled*pce/(pc_fit_scaled(end-1)*1000);
    krg_fit_scaled = krg_fit_scaled*krg_swc/krg_fit_scaled(1);
    
    % Assuming reversible krw, we have 3 points: at swc, at 1-sgt, and at
    % sw=1, so we fit a curve to find the exponent.
    swpts = [swc(n), 1-sgt, 1];
    nswpts = (swpts-swc(n)) ./ (1-swc(n));
    krpts = [krw_fit_scaled(1), krw_sgt, 1];
    w_fit = fittype('x^b', 'independent', 'x');
    w_model = fit(nswpts', krpts', w_fit);
%     if w_model.b < 1
%        disp('power krw exponent < 1, using linear model')
%        w_fit = fittype('a*x^1', 'independent', 'x');
%        w_model = fit(nswpts', krpts', w_fit);
%     end
    krw_fit_scaled = w_model(ns);
elseif model_case == 32
    pc_fit_scaled(:,n) = pc_fit_scaled1*pce(n)/(pc_fit_scaled1(end-1)*1000);
end

end

% Figure
cmap = copper(m);
h = figure(1);
latx  = {'Interpreter','latex'};
if plot_data
    subplot(1,5,2)
else
    subplot(1,3,1) 
end
hold on
l1 = plot(sw(:,1), pc_fit_scaled(:,1)*1000, '-', 'markersize', 4, ...
         'DisplayName', 'ESF', 'color', cmap(1,:), 'linewidth', 1);
l2 = plot(sw(:,2), pc_fit_scaled(:,2)*1000, '-', 'markersize', 4, ...
         'DisplayName', 'C', 'color', cmap(2,:), 'linewidth', 1);
l3 = plot(sw(:,3), pc_fit_scaled(:,3)*1000, '-', 'markersize', 4, ...
         'DisplayName', 'D', 'color', cmap(3,:), 'linewidth', 1);
l4 = plot(sw(:,4), pc_fit_scaled(:,4)*1000, '-', 'markersize', 4, ...
         'DisplayName', 'E', 'color', cmap(4,:), 'linewidth', 1);
if model_case < 3
    l5 = plot(sw(:,5), pc_fit_scaled(:,5)*1000, '-', 'markersize', 4, ...
             'DisplayName', 'F', 'color', cmap(5,:), 'linewidth', 1);
    leg = legend([l1 l2 l3 l4 l5], 'location', 'northeast', latx{:});
else
    leg = legend([l1 l2 l3 l4], 'location', 'northeast', latx{:});
end
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$p_\mathrm{c}$ [mbar]', latx{:}, 'fontSize', 14)
set(gca, 'yscale', 'log')
ylimv = 200;
ylim([0.1 ylimv]), yticks([0.1 0.5 1 5 10 50 100 ylimv])
xlim([0 1]), xticks(0:.2:1)
grid on
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.8]));
hold off
%title('p_c curves', 'fontweight', 'normal');

if plot_data
    subplot(1,5,4)
else
    subplot(1,3,2) 
end
if model_case < 3
    krg_name = '$k_{r\mathrm{g}}$: $0.97(1-S_w^*)^{1.39}$, Swc';
    krw_name = '$k_{r\mathrm{w}}$: $(S_w^*)^{4.192}$, Swc';
else
    krg_name = '$k_{r\mathrm{g}}$: $0.97(1-S_w^*)^{1.39}$, Swc, $k_{r\mathrm{g, swc}}$';
    krw_name = ['$k_{r\mathrm{w}}$: $(S_w^*)^b$, Swc'];
end
hold on
l1 = plot(sw(:,1), krg_fit_scaled(:,1), '-', 'linewidth', 1, ...
          'DisplayName', 'ESF', 'color', cmap(1,:));
l2 = plot(sw(:,2), krg_fit_scaled(:,2), '-', 'linewidth', 1, ...
          'DisplayName', 'C', 'color', cmap(2,:));
l3 = plot(sw(:,3), krg_fit_scaled(:,3), '-', 'linewidth', 1, ...
         'DisplayName', 'D', 'color', cmap(3,:));
if model_case < 3
    l4 = plot(sw(:,5), krg_fit_scaled(:,5), '-', 'linewidth', 1, ...
         'DisplayName', 'E, F, G', 'color', cmap(5,:));
else
    l4 = plot(sw(:,4), krg_fit_scaled(:,4), '-', 'linewidth', 1, ...
        'DisplayName', 'E', 'color', cmap(4,:));
    l5 = plot(sw(:,5), krg_fit_scaled(:,5), '-', 'linewidth', 1, ...
        'DisplayName', 'F', 'color', cmap(5,:));
    l6 = plot(sw(:,6), krg_fit_scaled(:,6), '-', 'linewidth', 1, ...
         'DisplayName', 'G', 'color', cmap(6,:));
end
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$k_\mathrm{rg}$ [-]', latx{:}, 'fontSize', 14)
%set(gca, 'yscale', 'log')
ylim([0 1]), yticks(0:.1:1)
xlim([0 1]), xticks(0:.2:1)
grid on
hold off
%title('k_{rg} curves', 'fontweight', 'normal');

if plot_data
    subplot(1,5,5)
else
    subplot(1,3,3) 
end
hold on
l1 = plot(sw(:,1), krw_fit_scaled(:,1), '-', 'linewidth', 1, ...
          'DisplayName', 'ESF', 'color', cmap(1,:));
l2 = plot(sw(:,2), krw_fit_scaled(:,2), '-', 'linewidth', 1, ...
          'DisplayName', 'C', 'color', cmap(2,:));
l3 = plot(sw(:,3), krw_fit_scaled(:,3), '-', 'linewidth', 1, ...
         'DisplayName', 'D', 'color', cmap(3,:));
if model_case < 3
    l4 = plot(sw(:,4), krw_fit_scaled(:,4), '-', 'linewidth', 1, ...
             'DisplayName', 'E, F, G', 'color', cmap(5,:));
else
    l4 = plot(sw(:,4), krw_fit_scaled(:,4), '-', 'linewidth', 1, ...
             'DisplayName', 'E', 'color', cmap(4,:));
    l5 = plot(sw(:,5), krw_fit_scaled(:,5), '-', 'linewidth', 1, ...
             'DisplayName', 'F', 'color', cmap(5,:));
    l6 = plot(sw(:,6), krw_fit_scaled(:,6), '-', 'linewidth', 1, ...
             'DisplayName', 'G', 'color', cmap(6,:));
    
end
hold off
xlabel('$S_\mathrm{w}$ [-]', latx{:}, 'fontSize', 14)
ylabel('$k_\mathrm{rw}$ [-]', latx{:}, 'fontSize', 14)
ylim([0 1]), yticks(0:.1:1)
xlim([0 1]), xticks(0:.2:1)
grid on
%title('k_{rw} curves', 'fontweight', 'normal');
if plot_data
    set(h, 'Position', [200, 200, 1500, 325])
else
    set(h, 'Position', [200, 200, 900, 325])
end


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
l11 = plot(data.sgi, data.sgt, 'sk', 'DisplayName', 'Pentland et al. (2010)');
l22 = plot(sgi, sgt_vals, '--k', 'DisplayName', ...
          '$S_\mathrm{gt}^* = S_\mathrm{gi}^*/(1+5.2S_\mathrm{gi}^*)$');
grid on
hold off
xlabel('S_{gi}')
ylabel('S_{gt}')
xlim([0 1]), xticks(0:.2:1)
ylim([0 0.15]), yticks(0:0.025:0.15)
legend([l11, l22(1)], 'location', 'southeast', latx{:})
set(h, 'Position', [200, 200, 350, 275])

% Calculate imbibition curves
% Gas (Land, 1968)
sgt_max = zeros(1,m);
for n=1:m
    sgi_max = max(sgi(:,n));
    nsgi_max = sgi_max / (1-swc(n)); % this is always 1 for bounding ib curve
    if model_case < 3
        nsgt_max = sgt_model(nsgi_max);
        sgt_max(n) = nsgt_max*(1-swc(n));
    elseif model_case == 3 || model_case == 32
        sgt_max(n) = sgt(n);
        nsgt_max = sgt_max(n)/(1-swc(n));
    end
    %nsgt_max = 0.4;                              % Land (1968) Fig. 2
    C(n) = 1/nsgt_max - 1;                           % Land (1968) Eq. 2
    sgib(:,n) = linspace(sgt_max(n), sgi_max, np);
    nsgib = (sgib(:,n))/(1-swc(n));            
    nsgibF = 0.5*((nsgib-nsgib(1)) + sqrt((nsgib-nsgib(1)).^2 + ...
                                          (4/C(n))*(nsgib-nsgib(1))));  % Land (1968) Eq. 4
    e = 2/bc_model.lambda + 3;
    %e = 4;                                      % Land (1968) Fig. 2
    krgib1 = nsgibF.^2.*(1-(1-nsgibF).^(e-2));
    krgib(:,n) = max(krg_fit_scaled(:,n))*krgib1;             % our krg is not 1 at sgi;

    % Water and Pcog (no hysteresis)
    nswib = 1 - nsgib;
    krwib(:,n) = w_model(nswib);
    nspcfitib = [nswib(1:end-1); nspcmin];
    pc_fit_scaled_ib = interp1(nspcfit, pc_fit_scaled, nspcfitib)';
end

figure(1)
if plot_data
    subplot(1,5,4)
else
    subplot(1,3,2) 
end
hold on
l7 = plot(1-sgib(:,1), krgib(:,1), '--', 'linewidth', 1, ...
          'DisplayName', '$k_{r\mathrm{g}}^\mathrm{ib}$', 'color', cmap(1,:));
if model_case < 3
    plot(1-sgib(:,6), krgib(:,6), '--', 'linewidth', 1, ...
     'color', cmap(5,:))
    leg = legend([l1, l2, l3, l4, l7], 'location', 'northeast', latx{:});
else
    plot(1-sgib(:,2), krgib(:,2), '--', 'linewidth', 1, ...
     'color', cmap(2,:))
    plot(1-sgib(:,3), krgib(:,3), '--', 'linewidth', 1, ...
     'color', cmap(3,:))
    plot(1-sgib(:,4), krgib(:,4), '--', 'linewidth', 1, ...
     'color', cmap(4,:))
    plot(1-sgib(:,5), krgib(:,5), '--', 'linewidth', 1, ...
     'color', cmap(5,:))
    plot(1-sgib(:,6), krgib(:,6), '--', 'linewidth', 1, ...
     'color', cmap(6,:))
    leg = legend([l1, l2, l3, l4, l5, l6, l7], 'location', 'northeast', latx{:});
end
hold off
set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.8]));
set(gcf,'PaperType','A2')
%figure(1); subplot(1,3,2); xticks([0.84:.02:0.94]); xlim([0.84 0.94])
%yticks([0:.002:.008]); ylim([0 0.008])
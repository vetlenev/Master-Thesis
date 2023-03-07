function comparisonPlot(faults, SSF)
%
%
%

% Get kxx 
kxx = cellfun(@(x) x.Perm(1), faults)/(milli*darcy);
K = log10(kxx);
nbins = 25;
nedg = nbins + 1;
logMinP = min(min(K));
logMaxP = max(max(K));
edges = linspace(fix(logMinP)-1, fix(logMaxP)+1, nedg);

% Plot histograms
latx = {'Interpreter', 'latex'};
sz = [13, 12];
fh = figure(randi(90000, 1, 1));
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
labl = "$\log_{10}(k_{xx}$ [mD])";
nexttile(1)
histogram(K(:, 4), edges, 'Normalization', 'probability', ...
         'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlabel(labl(1), latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP) fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
xticks(fix(logMinP):1:fix(logMaxP)+1)
text(1, 0.5, ['SSF = ' num2str(SSF(4))], latx{:}, 'fontsize', sz(1))

nexttile(2)
histogram(K(:, 3), edges, 'Normalization', 'probability', ...
         'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlim([fix(logMinP) fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
xticks(fix(logMinP):1:fix(logMaxP)+1)
text(1, 0.5, ['SSF = ' num2str(SSF(3))], latx{:}, 'fontsize', sz(1))

nexttile(3)
histogram(K(:, 2), edges, 'Normalization', 'probability', ...
         'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlim([fix(logMinP) fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
xticks(fix(logMinP):1:fix(logMaxP)+1)
text(1, 0.5, ['SSF = ' num2str(SSF(2))], latx{:}, 'fontsize', sz(1))

nexttile(4)
histogram(K(:, 1), edges, 'Normalization', 'probability', ...
         'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlim([fix(logMinP) fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
xticks(fix(logMinP):1:fix(logMaxP)+1)
text(1, 0.5, ['SSF = ' num2str(SSF(1))], latx{:}, 'fontsize', sz(1))
set(fh, 'position', [200, 200, 175, 500]);

fh = figure(randi(90000, 1, 1));
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
labl = "$\log_{10}(k_{xx}$ [mD])";
nexttile(1)
histogram(K(:, 7), edges, 'Normalization', 'probability', ...
         'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlabel(labl(1), latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP) fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
xticks(fix(logMinP):1:fix(logMaxP)+1)
text(1, 0.5, ['SSF = ' num2str(SSF(7))], latx{:}, 'fontsize', sz(1))

nexttile(2)
histogram(K(:, 6), edges, 'Normalization', 'probability', ...
         'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlim([fix(logMinP) fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
xticks(fix(logMinP):1:fix(logMaxP)+1)
text(1, 0.5, ['SSF = ' num2str(SSF(6))], latx{:}, 'fontsize', sz(1))

nexttile(3)
histogram(K(:, 5), edges, 'Normalization', 'probability', ...
         'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlim([fix(logMinP) fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
grid on
xticks(fix(logMinP):1:fix(logMaxP)+1)
text(1, 0.5, ['SSF = ' num2str(SSF(5))], latx{:}, 'fontsize', sz(1))
set(fh, 'position', [400, 200, 175, 375]);

% Compute permeabilities from Kettermann et al., JGR:SE Fig. 14 (Exp 3 and
% 8, which is more similar to shear failure)
to_m2s = 1/(10^4*60);                                       % cm2/min to m2/s
exp3 = [1 1.5 2:10;                                         % SSF
        0.23 0.5 0.9 1.7 3.4 3.9 4.9 5.75 6.5 7.6 7.53];    % Q/dh [cm2/min]
exp3(2,:) = exp3(2,:)*to_m2s;
exp4 = [1:0.5:4 5:10;
        0.16 3 5.1 6.6 7.4 8.1 8.8 9.9 10.7 11.2 11.7 12.1 12.60];
exp4(2,:) = exp4(2,:)*to_m2s;
exp7 = [1:0.5:3 4:10;
        0.9 2 2.6 3.2 3.9 5.1 5.7 6.2 7 8.1 7.9 9.03];
exp7(2,:) = exp7(2,:)*to_m2s;        
exp8 = [1:0.5:4 5:10;
        0.21 0.3 0.45 0.7 1.1 1.6 2.2 3.3 3.8 4.4 5.3 6.2 7.03];
exp8(2,:) = exp8(2,:)*to_m2s;  

% Darcy's law params
A = 0.5*0.3;                        % outflow base area (m2)
rho = 998.23;                       % water density at 20C [kg/m3]
g = 9.81;                           % gravity
mu = 10^-3;                         % water viscosity at 20C [Pa*s]
% L3 = 0.28 - exp3(1,:)/100;          % Length of porous media (m), acciounting for throw
% L4 = 0.28 - exp4(1,:)/100;
% L7 = 0.28 - exp7(1,:)/100;
% L8 = 0.28 - exp8(1,:)/100;
[L3, L4, L7, L8] = deal(0.28);

% Hydraulic conductivity [L/T]
K3 = exp3(2,:).*L3/A;
K4 = exp4(2,:).*L4/A;
K7 = exp7(2,:).*L7/A;
K8 = exp8(2,:).*L8/A;

% Permeability [L2]
k3 = K3*mu/(rho*g)/(milli*darcy);
k4 = K4*mu/(rho*g)/(milli*darcy);
k7 = K7*mu/(rho*g)/(milli*darcy);
k8 = K8*mu/(rho*g)/(milli*darcy);

% Plot
ks = [0.6 1 1.4]*10^-4*mu/(rho*g)/(milli*darcy);
kc = [4 2 1]*10^-8*mu/(rho*g)/(milli*darcy);
cmap = copper(10);
latx = {'Interpreter', 'latex'};
figure(39339)
subplot(1,2,1)
hold on
p1 = fill(log10([ks(1) ks(end) ks(end) ks(1) ks(1)]), [1 1 10 10 1], ...
          cmap(end-2, :), 'faceAlpha', 0.3, 'edgecolor', 'none', ...
          'DisplayName', '$\pm e$');
fill(log10([kc(1) kc(end) kc(end) kc(1) kc(1)]), [1 1 10 10 1], ...
     cmap(3, :), 'faceAlpha', 0.3, 'edgecolor', 'none')
p2 = plot(log10([ks(2) ks(2)]), [1 10], '-', 'color', cmap(end-2, :), ...
          'linewidth', 1, 'DisplayName', 'sand');
p3 = plot(log10([kc(2) kc(2)]), [1 10], '-', 'color', cmap(3, :), ...
          'linewidth', 1, 'DisplayName', 'clay');
p4 = plot(log10(k3), exp3(1,:), 'ok', 'markerFaceColor', [0.2 0.2 0.2], ...
          'DisplayName', 'exp 3');
p5 = plot(log10(k4), exp4(1,:), 'sk', 'markerFaceColor', [0.45 0.45 0.45], ...
          'DisplayName', 'exp 4');
ylim([1 10])
xlim([0 4.5])
xlabel('$\log_{10}(k$ [mD])', 'fontSize', 12, latx{:})
ylabel('$t/T$ (SSF) [-]', 'fontSize', 12, latx{:})
title('Single clay layer', 'fontsize', 14, latx{:})
grid on
h = legend([p1 p2 p3 p4 p5], latx{:}, 'location', 'northwest', 'fontsize', 10);
set(h.BoxFace, 'ColorType','truecoloralpha', ...
    'ColorData', uint8(255*[1;1;1;.5]));

subplot(1,2,2)
hold on
fill(log10([ks(1) ks(end) ks(end) ks(1) ks(1)]), [1 1 10 10 1], ...
          cmap(end-2, :), 'faceAlpha', 0.3, 'edgecolor', 'none', ...
          'DisplayName', '$\pm e$')
fill(log10([kc(1) kc(end) kc(end) kc(1) kc(1)]), [1 1 10 10 1], ...
     cmap(3, :), 'faceAlpha', 0.3, 'edgecolor', 'none')
plot(log10([ks(2) ks(2)]), [1 10], '-', 'color', cmap(end-2, :), ...
     'linewidth', 1, 'DisplayName', 'sand')
plot(log10([kc(2) kc(2)]), [1 10], '-', 'color', cmap(3, :), ...
     'linewidth', 1, 'DisplayName', 'clay')
p4 = plot(log10(k7), exp7(1,:), 'dk', 'markerFaceColor', [0.65 0.65 0.65], ...
          'DisplayName', 'exp 7');
p5 = plot(log10(k8), exp8(1,:), 'pk', 'markerFaceColor', [0.9 0.9 0.9], ...
          'DisplayName', 'exp 8');
ylim([1 10])
xlim([0 4.5])
xlabel('$\log_{10}(k$ [mD])', 'fontSize', 12, latx{:})
ylabel('$t/T$ (SSF) [-]', 'fontSize', 12, latx{:})
title('Double clay layer', 'fontsize', 14, latx{:})
grid on
h = legend([p4 p5], latx{:}, 'location', 'northwest', 'fontsize', 10);
set(h.BoxFace, 'ColorType','truecoloralpha', ...
    'ColorData', uint8(255*[1;1;1;.5]));


end
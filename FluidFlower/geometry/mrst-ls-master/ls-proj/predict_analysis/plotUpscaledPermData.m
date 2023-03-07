function plotUpscaledPermData(faults, thickness, vcl, zf, zmax, name, dim)
%
%
%

% Utilities
latx = {'Interpreter', 'latex'};
sz = [14, 12];

% Fault MatProps
if dim == 2
    kAlongStrike = cell2mat(cellfun(@(x) x.Grid.permy, faults, ...
                                    'UniformOutput', false)')./(milli*darcy);
end
perms = cell2mat(cellfun(@(x) x.Perm, faults, ...
                         'UniformOutput', false)) ./ (milli*darcy);
if dim == 2
    thick = cell2mat(cellfun(@(x) x.MatProps.thick, faults, ...
                             'UniformOutput', false));
elseif dim == 3
    thick = cell2mat(cellfun(@(x) x.Thick, faults, 'UniformOutput', false));
end
fvcl = cell2mat(cellfun(@(x) x.Vcl, faults, 'UniformOutput', false));
if any(any(perms < 0))
    id = unique([find(perms(:, 1)<0), find(perms(:, 2)<0), find(perms(:, 3)<0)]);
    warning(['Negative upscaled perms found in ' num2str(numel(id))...
             ' simulations (ignored).'])
    kAlongStrike(:, id) = [];
    perms(id, :) = [];
    thick(id, :) = [];
    fvcl(id, :) = [];
end
if dim == 2
    logkStrikeBounds = log10(getAveragingPerm(kAlongStrike, {'ha', 'ah'}));
end
    
% Hist params
K = log10(perms);
nedg = 26;
logMinP = min(min(K));
logMaxP = max(max(K));
edgesk = linspace(fix(logMinP)-1, fix(logMaxP)+1, nedg);

% SGR and spe02 perm
edges = linspace(0.1, 0.6, nedg)';
SGR = getSGR(thickness, vcl);
SGR = [min(SGR(:, 1)), mean(SGR(:, 1)), max(SGR(:, 1))];
logk_spe02 = log10(getPermSpe02(SGR, zf(1), zmax{1}(1)));

% Literature measurements
if strcmp(name, 'GoM')
    % GoM data (GoM Atlas, 2017, Ch. 3)
    % XRD (Mass fraction assumed \approx Vcl since rho_s_c and rho_s_bulk are 
    % probably similar)
    z = [10578 10609]*0.3048;  
    mcl = [0.27 0.38];
    n = [0.0315 0.0932];
    logk_lit = log10([0.000145 0.00455]);      % log(mD)
    cf_lit = mcl;
    dispName = 'Lu17';
    linSty = '--';
elseif strcmp(name, 'Tar')
    % Taranaki data (Childs et al., 2007)
    vclt = [0.225 0.51];
    logk_lit = [-3.699 0];
    cf_lit = vclt;
    dispName = 'Chi07';
    linSty = '-.';
end

%% Perm Histograms    
fh = figure(randi(10000, 1, 1));
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
labls = ["$\log_{10}(k_{xx}$ [mD])", ...
    "$\log_{10}(k_{yy}$ [mD])", ...
    "$\log_{10}(k_{zz}$ [mD])"];

nexttile(1)
hold on
l1 = fill([logk_spe02(3) logk_spe02(1) logk_spe02(1) logk_spe02(3) logk_spe02(3)], ...
     [0.001 0.001 0.6 0.6 0.6], ...
     [0.85 0.85 0.85], 'edgecolor', 'none', 'DisplayName', 'Spe02');
if strcmp(name, 'GoM') || strcmp(name, 'Tar')
    l2 = plot([logk_lit(1) logk_lit(1)], [0.001 0.6], ...
              'k', 'lineStyle', linSty, 'lineWidth', 1.25, 'DisplayName', dispName);
    plot([logk_lit(2) logk_lit(2)], [0.001 0.6], ...
         'k', 'lineStyle', linSty, 'lineWidth', 1.25);
end
histogram(K(:, 1), edgesk, 'Normalization', 'probability', ...
          'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
xlabel(labls(1), latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:.5)
grid on
% if strcmp(name, 'GoM') || strcmp(name, 'Tar')
%     leg = legend([l1 l2], latx{:}, 'fontSize', 9);
% else
%     leg = legend(l1, latx{:}, 'fontSize', 9);
% end
% set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))

nexttile(2)
rr = [255, 125, 125]/255;
if dim == 2
    plot(repelem(logkStrikeBounds(1), 2), [0 1], '-', ...
        'color', 'k', 'lineWidth', 1);
    hold on
    plot(logkStrikeBounds(1), 0.25, 'dk', 'markerFacecolor', rr, 'markerSize', 4)
    plot(repelem(logkStrikeBounds(2), 2), [0 1], '-', ...
        'color', 'k', 'lineWidth', 1);
    plot(logkStrikeBounds(2), 0.25, 'dk', 'markerFacecolor', rr, 'markerSize', 4)
end
histogram(K(:, 2), edgesk, 'Normalization', 'probability', ...
    'FaceColor', rr, 'FaceAlpha', 1)
xlabel(labls(2), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:.5)
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
hold off

nexttile(3)
bb = [125, 125, 255]/255;
histogram(K(:, 3), edgesk, 'Normalization', 'probability', ...
    'FaceColor', bb, 'FaceAlpha', 1)
xlabel(labls(3), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:.5)
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
set(fh, 'position', [200, 200, 150, 350]);


%% Correlations and clay Fraction histogram
fh = figure(randi(10000, 1, 1));
if dim == 2
    % fault thickness  vs perm
    [Rx, Px] = corr(K(:,1), thick, 'Type', 'Spearman');         % corrcoeff and pval matrices
    [Ry, Py] = corr(K(:,2), thick, 'Type', 'Spearman');
    [Rz, Pz] = corr(K(:,3), thick, 'Type', 'Spearman');
    a = 0.05;                                    % significance level
    pvals = [Px Py Pz];
    r = [Rx Ry Rz];
    if pvals(1) < a, colr{1} = 'k'; fw{1} = 'bold';
    else, colr{1} = 'k'; fw{1} = 'normal'; end
    if pvals(2) < a, colr{2} = 'r'; fw{2} = 'bold';
    else, colr{2} = 'r'; fw{2} = 'normal'; end
    if pvals(3) < a, colr{3} = 'b'; fw{3} = 'bold';
    else, colr{3} = 'b'; fw{3} = 'normal'; end

    tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    nexttile(3)
    hold on
    scatter(thick, K(:,1), 8, 'ok','MarkerEdgeAlpha', 0.6)
    scatter(thick, K(:,2), 8, 'sr','MarkerEdgeAlpha', 0.3)
    scatter(thick, K(:,3), 8, 'db','MarkerEdgeAlpha', 0.15)
    if any(pvals ~= 0) && any(log10(pvals) > -10)
        title({['$r_s = $ ' num2str(round(r(1), 2)) ', ' ...
            num2str(round(r(2), 2)) ', ' num2str(round(r(3), 2))]; ...
            ['$\log_{10} (p) = $' num2str(round(log10(pvals(1)))) ', ' ...
            num2str(round(log10(pvals(2)))) ', ' ...
            num2str(round(log10(pvals(3))))]}, ...
            latx{:}, 'fontSize', 10);
    else
        title({['$r_s = $ ' num2str(round(r(1), 2)) ', ' ...
            num2str(round(r(2), 2)) ', ' num2str(round(r(3), 2))]}, ...
            latx{:}, 'fontSize', 10);
    end
    xlabel('f$_\mathrm{T}$ [m]', latx{:}, 'fontSize', sz(2))
    ylabel('$\log_{10}(k$ [mD])', latx{:}, 'fontSize', sz(2))
    ylim([fix(logMinP)-1 fix(logMaxP)+1])
    xlim([0 fix(max(thick))+1])
    h = legend({'$k_{xx}$', '$k_{yy}$', '$k_{zz}$'}, latx{:}, 'fontSize', 10);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
        'ColorData', uint8(255*[1;1;1;.7]));
    grid on
    hold off
else
    tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
end

% N smear domains vs perm
if dim == 2
    type = 'Spearman';
else
    type = 'Pearson';
end
[Rx, Px] = corr(K(:,1), fvcl, 'Type', type);         % corrcoeff and pval matrices
[Ry, Py] = corr(K(:,2), fvcl, 'Type', type);
[Rz, Pz] = corr(K(:,3), fvcl, 'Type', type);
a = 0.05;                                    % significance level
pvals = [Px Py Pz];
r = [Rx Ry Rz];
if pvals(1) < a, colr{1} = 'm'; fw{1} = 'bold';
else, colr{1} = 'k'; fw{1} = 'normal'; end
if pvals(2) < a, colr{2} = 'm';
else, colr{2} = 'r'; fw{2} = 'normal'; end
if pvals(3) < a, colr{3} = 'm';
else, colr{3} = 'b'; fw{3} = 'normal'; end
nexttile(2)
hold on
scatter(fvcl, K(:,1), 8, 'ok','MarkerEdgeAlpha', 0.6)
scatter(fvcl, K(:,2), 8, 'sr','MarkerEdgeAlpha', 0.3)
scatter(fvcl, K(:,3), 8, 'db','MarkerEdgeAlpha', 0.15)
if dim == 2
    corrtype = '$r_s = $ ';
else
    corrtype = '$\rho = $ ';
end
if any(pvals ~= 0) && any(log10(pvals) > -10)
    title({[corrtype num2str(round(r(1), 2)) ', ' ...
        num2str(round(r(2), 2)) ', ' num2str(round(r(3), 2))]; ...
        ['$\log_{10} (p) = $' num2str(round(log10(pvals(1)))) ', ' ...
        num2str(round(log10(pvals(2)))) ', ' ...
        num2str(round(log10(pvals(3))))]}, ...
        latx{:}, 'fontSize', 10);
else
    title({[corrtype num2str(round(r(1), 2)) ', ' ...
        num2str(round(r(2), 2)) ', ' num2str(round(r(3), 2))]}, ...
        latx{:}, 'fontSize', 10);
end
xlabel('f$_{V_\mathrm{cl}}$ [-]', latx{:}, 'fontSize', sz(2))
ylabel('$\log_{10}(k$ [mD])', latx{:}, 'fontSize', sz(2))
ylim([fix(logMinP)-1 fix(logMaxP)+1])
xlim([min(fvcl) max(fvcl)])
%xticks([0 0.2 0.4 0.6 0.8 1])
if dim == 3
    h = legend({'$k_{xx}$', '$k_{yy}$', '$k_{zz}$'}, latx{:}, 'fontSize', 10);
    set(h.BoxFace, 'ColorType','truecoloralpha', ...
        'ColorData', uint8(255*[1;1;1;.7]));
end
grid on
hold off

% SGR histogram
clrs = copper(10);
if strcmp(name, 'GoM')
    idclr = 3;
elseif strcmp(name, 'Tar')
    idclr = 7;
else
    idclr = 5;
end
nexttile(1)
hold on
l1 = fill([SGR(1) SGR(3) SGR(3) SGR(1) SGR(1)], [0.001 0.001 0.499 0.499 0.001], ...
     [255, 230, 196]/255, 'edgecolor', 'none', 'FaceAlpha', 0.8, 'DisplayName', 'SGR');
if strcmp(name, 'GoM') || strcmp(name, 'Tar')
    l3 = plot([cf_lit(1) cf_lit(1)],[0.001 0.5], 'k', 'lineStyle', linSty, 'lineWidth', 1.25, ...
              'DisplayName', 'Lu17');
    plot([cf_lit(2) cf_lit(2)],[0.001 0.5], 'k', 'lineStyle', linSty, 'lineWidth', 1.25);
end
histogram(fvcl(:, 1), edges, 'Normalization', 'probability', ...
          'FaceColor', clrs(idclr, :), 'FaceAlpha', 1)
xlabel('f$_{V_\mathrm{cl}}$ [-]', latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xticks([0.1 0.2 0.3 0.4 0.5 0.6]); xlim([0.1 0.6])
ylim([0 0.5]); yticks(0:.1:.5)
grid on
%leg = legend([l1 l2 l3], latx{:}, 'fontSize', 9);
%set(leg.BoxFace, 'ColorType','truecoloralpha', 'ColorData', uint8(255*[1;1;1;.5]));
if dim == 2
    set(fh, 'position', [200, 200, 175, 400]);
else
    set(fh, 'position', [200, 200, 150, 250]);
end

end
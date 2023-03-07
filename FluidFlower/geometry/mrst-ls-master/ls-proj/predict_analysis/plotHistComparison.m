function plotHistComparison(faults, thickness, vcl, zf, zmax)
%
%
%

% Utilities
latx = {'Interpreter', 'latex'};
sz = [18, 18, 14];

% Fault MatProps
kAlongStrike = cell2mat(cellfun(@(x) x.Grid.permy, faults, ...
                                'UniformOutput', false)')./(milli*darcy);
perms = cell2mat(cellfun(@(x) x.Perm, faults, ...
                         'UniformOutput', false)) ./ (milli*darcy);
thick = cell2mat(cellfun(@(x) x.MatProps.thick, faults, ...
                         'UniformOutput', false));
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
logkStrikeBounds = log10(getAveragingPerm(kAlongStrike, {'ha', 'ah'}));
    
% Hist params
K = log10(perms);
nedg = 26;
logMinP = min(min(K));
logMaxP = max(max(K));
edges = linspace(fix(logMinP)-1, fix(logMaxP)+1, nedg);

% Spe02 calculation
SGR = mean(getSGR(thickness, vcl));
logSpe02 = log10(getPermSpe02(SGR, zf(1,1), zmax{1}(1)));
if zf(1,1) == 500 % shallow fault
    krat = 3;
else              % medium depth fault
    krat = 10;
end

% LDL calculation
A = log10(getPermSpe02(0, zf(1,1), zmax{1}(1)));
if zf(1,1) == 500
    B = -0.15; c = 0; d = 5; E = -0.02;
else
    B = -0.1; c = 0; d = 3.5; E = -0.012;
end

SGR = 100*SGR;
if SGR < 15
    logGra19 = A - c + SGR*(B*15-(d-c))/15;
elseif SGR >= 15 && SGR <= 40
    logGra19 = A - d + B*SGR;
else
    logGra19 = A - d + 40*B + E*(SGR-40);
end


% --------------------------- Histograms ----------------------------------
fh = figure(randi(10000, 1, 1));
tiledlayout(1, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
labls = ["$k_{xx}$", "$k_{yy}$", "$k_{zz}$"];
lablx = "$\log_{10}(k$ [mD])";
nexttile(1)
hold on
plot(repelem(min(K(:,1)), 2), [0 1], '--', 'color', [0.5 0.5 0.5], 'lineWidth', 2)
plot(repelem(max(K(:,1)), 2), [0 1], '--', 'color', [0.5 0.5 0.5], 'lineWidth', 2)
plot(repelem(logSpe02, 2), [0 1], '-', 'color', 'k', 'lineWidth', 1);
l1 = plot(logSpe02, 0.5, 'ok', 'markerFacecolor', [0.8 0.8 0.8], ...
          'markerSize', 8, 'DisplayName', 'Spe02');
plot(repelem(logGra19, 2), [0 1], '-', 'color', 'k', 'lineWidth', 1);
l2 = plot(logGra19, 0.5, 'sk', 'markerFacecolor', 'y', 'markerSize', 9, ...
          'DisplayName', 'LDL');
histogram(K(:, 1), edges, 'Normalization', 'probability', ...
          'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1)
title(labls(1), 'fontSize', sz(2), latx{:})
xlabel(lablx, latx{:}, 'fontSize', sz(2))
ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
h = gca; h.FontSize = sz(3); 
xticks(fix(logMinP)-1:2:fix(logMaxP)+1);
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
%leg = legend([l1(1), l2(1)], 'location', 'northwest');
%set(leg.BoxFace, 'ColorType','truecoloralpha', ...
%    'ColorData', uint8(255*[1;1;1;.5]));

nexttile(2)
rr = [255, 125, 125]/255;
hold on
plot(repelem(min(K(:,2)), 2), [0 1], '--', 'color', rr, 'lineWidth', 2)
plot(repelem(max(K(:,2)), 2), [0 1], '--', 'color', rr, 'lineWidth', 2)
plot(repelem(logkStrikeBounds(1), 2), [0 1], '-', ...
    'color', 'k', 'lineWidth', 1);
plot(logkStrikeBounds(1), 0.5, 'dk', 'markerFacecolor', rr, 'markerSize', 6)
plot(repelem(logkStrikeBounds(2), 2), [0 1], '-', ...
    'color', 'k', 'lineWidth', 1);
plot(logkStrikeBounds(2), 0.5, 'dk', 'markerFacecolor', rr, 'markerSize', 6)
histogram(K(:, 2), edges, 'Normalization', 'probability', ...
    'FaceColor', rr, 'FaceAlpha', 1)
title(labls(2), 'fontSize', sz(2), latx{:})
%xlabel(labls(2), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
xticks(fix(logMinP)-1:2:fix(logMaxP)+1);
set(gca, 'xticklabel', {[]}, 'yticklabel', {[]})
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
hold off

nexttile(3)
bb = [125, 125, 255]/255;
hold on
plot(repelem(min(K(:,3)), 2), [0 1], '--', 'color', bb, 'lineWidth', 2)
plot(repelem(max(K(:,3)), 2), [0 1], '--', 'color', bb, 'lineWidth', 2)
plot(repelem(log10((10^logSpe02)*krat), 2), [0 1], '-', ...
     'color', 'k', 'lineWidth', 1);
plot(log10((10^logSpe02)*krat), 0.5, 'ok', 'markerFacecolor', [0.8 0.8 0.8], 'markerSize', 8)
plot(repelem(log10((10^logGra19)*krat), 2), [0 1], '-', 'color', 'k', 'lineWidth', 1);
plot(log10((10^logGra19)*krat), 0.5, 'sk', 'markerFacecolor', 'y', 'markerSize', 9)
histogram(K(:, 3), edges, 'Normalization', 'probability', ...
    'FaceColor', bb, 'FaceAlpha', 1)
title(labls(3), 'fontSize', sz(2), latx{:})
%xlabel(labls(3), latx{:}, 'fontSize', sz(2))
%ylabel('P [-]', latx{:}, 'fontSize', sz(2))
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.6]); yticks(0:.2:.6)
xticks(fix(logMinP)-1:2:fix(logMaxP)+1);
set(gca, 'xticklabel', {[]}, 'yticklabel', {[]})
grid on
%xticks(10.^(fix(logMinP)-1:2:fix(logMaxP)+1))
set(fh, 'position', [200, 200, 450, 175]);

end
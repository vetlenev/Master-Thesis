function plotMultipleBoxPlot(faults_all, Nsim, idSeq, thickness, vcl, ...
                             zf, zmax, opt)
%
%
%

% Find modes
nStrat = numel(faults_all);
Kxx = zeros(Nsim, nStrat);
Kzz = zeros(Nsim, nStrat);
logLim = zeros(nStrat, 2);
fVcl = zeros(Nsim, nStrat);
for n=1:nStrat
    K = log10(cell2mat(cellfun(@(x) x.Perm, faults_all{n}, ...
              'UniformOutput', false))./(milli*darcy));
    fVcl(:, n) = cell2mat(cellfun(@(x) x.Vcl, faults_all{n}, ...
                          'UniformOutput', false));
    Kxx(:, n) = K(:, 1);
    Kzz(:, n) = K(:, end);
    logLim(n, :) = [min(min(Kxx)), max(max(Kzz))];
end
nedg = 26;
edges = linspace(fix(min(logLim(:,1)))-1, fix(max(logLim(:,2)))+1, nedg);
locs = edges(1:end-1)+diff(edges)/2;
Pxx = zeros(nedg-1, nStrat);
Pzz = zeros(nedg-1, nStrat);
idModex = cell(nStrat, 1);
idModez = cell(nStrat, 1);
for n=1:nStrat
    Pxx(:, n) = histcounts(Kxx(:, n), edges, 'Normalization','probability');
    Pzz(:, n) = histcounts(Kzz(:, n), edges, 'Normalization','probability');
    idModex{n} = islocalmax(Pxx(:,n), 'MinSeparation', 1, ...
                            'SamplePoints', locs, 'MinProminence', 0.005);
    idModez{n} = islocalmax(Pzz(:,n), 'MinSeparation', 1, ...
                            'SamplePoints', locs, 'MinProminence', 0.005);
    if sum(idModex{n}) == 0
        [~, id] = max(Pxx(:, n));
        idModex{n}(id) = true;
    end
    if sum(idModez{n}) == 0
        [~, id] = max(Pzz(:, n));
        idModez{n}(id) = true;
    end
end

% check
% figure(683)
% ids = 14;
% subplot(1,2,1)
% histogram(Kxx(:, ids), edges, 'Normalization', 'probability')
% hold on
% plot(locs(idModex{ids}), Nxx(idModex{ids}, ids), 'or')
% subplot(1,2,2)
% histogram(Kzz(:, ids), edges, 'Normalization', 'probability')
% hold on
% plot(locs(idModez{ids}), Nzz(idModez{ids}, ids), 'or')

%% Boxplot plot 1
p595x = prctile(Kxx, [5 95])';      p595z = prctile(Kzz, [5 95])';
pMinMaxx = [min(Kxx); max(Kxx)]';   pMinMaxz = [min(Kzz); max(Kzz)]';
latx = {'Interpreter', 'latex'};
sz = [14, 12];
fh = figure(1);
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile(1)
l1 = plot(pMinMaxx, [1:nStrat; 1:nStrat]','*r','markersize', 3, ...
          'DisplayName','Min-Max');
hold on
boxplot(Kxx, idSeq, 'orientation','horizontal');
% Outliers
h2=findobj(gca,'tag','Outliers');
delete(h2)
h = flipud(findobj(gca,'Tag','Upper Whisker'));
set(h,'linewidth',1,'linestyle','-')
for k=1:length(h)
    ydata = get(h(k),'XData');
    ydata(2) = p595x(k,2);
    set(h(k),'XData',ydata);
end
h = flipud(findobj(gca,'Tag','Upper Adjacent Value'));
set(h,'linewidth',1,'marker','none')
%h.CapSize = 12
for k=1:length(h)
    ydata = get(h(k),'XData');
    ydata(:) = p595x(k,2);
    set(h(k),'XData',ydata);
end
h = flipud(findobj(gca,'Tag','Lower Whisker'));
set(h,'linewidth',1,'linestyle','-')
for k=1:length(h)
    ydata = get(h(k),'XData');
    ydata(1) = p595x(k, 1);
    set(h(k),'XData',ydata);
end
h = flipud(findobj(gca,'Tag','Lower Adjacent Value'));
set(h,'linewidth',1,'marker','none')
for k=1:length(h)
    ydata = get(h(k),'XData');
    ydata(:) = p595x(k, 1);
    set(h(k),'XData',ydata);
end

if strcmp(opt, 'full')
    for n=1:nStrat
        % Mode add
        [~, idMainMode] = max(Pxx(idModex{n}, n)); 
        modes = locs(idModex{n});
        l2 = plot(modes, n, 'sk', 'MarkerFaceColor', [1 1 0.7], ...
                  'MarkerSize', 7, 'DisplayName','Mode');
        
        l21 = plot(modes(idMainMode), n, 'sk', 'MarkerFaceColor', 'y', ...
                   'MarkerSize', 8, 'DisplayName','Main mode');
        
        % Spe02 add
        SGR = mean(getSGR(thickness{n}, vcl{n}));
        logSpe02 = log10(getPermSpe02(SGR, zf(n,1), zmax{n}{1}(1)));
        l3 = plot(logSpe02, n, 'ok', 'MarkerFaceColor', [0.8 0.8 0.8], ...
            'MarkerSize', 4, 'DisplayName','Spe02');
        
        % Grant19 LDL add
        A = log10(getPermSpe02(0, zf(n,1), zmax{n}{1}(1)));
        if mod(n, 2) == 0
            B = -0.15;
            c = 0;
            d = 5;
            E = -0.02;
        else
            B = -0.1;
            c = 0;
            d = 3.5;
            E = -0.012;
        end
        
        SGR = 100*SGR;
        if SGR < 15
            logGra19 = A - c + SGR*(B*15-(d-c))/15;
        elseif SGR >= 15 && SGR <= 40
            logGra19 = A - d + B*SGR;
        else
            logGra19 = A - d + 40*B + E*(SGR-40);
        end
        l4 = plot(logGra19, n, 'dk', 'MarkerFaceColor', [0, 138, 255]/255, ...
            'MarkerSize', 3, 'DisplayName','LDL');
    end
    ylabel('$\overline{V_\mathrm{cl}}$ [-]', 'fontSize', sz(2), latx{:})
    xlabel('$\log_{10}(k_{xx}$ [mD])', 'fontSize', sz(2), latx{:})
    grid on
    leg = legend([l1(1), l2(1), l21(1), l3, l4]);
    set(leg.BoxFace, 'ColorType','truecoloralpha', ...
        'ColorData', uint8(255*[1;1;1;.5]));
    
end


%% Boxplot plot 2
nexttile(2)
l1 = plot(pMinMaxz, [1:nStrat; 1:nStrat]','*r','markersize', 3, ...
          'DisplayName','Min-Max');
hold on
boxplot(Kzz, idSeq, 'orientation','horizontal');
% Outliers
h2=findobj(gca,'tag','Outliers');
delete(h2)
h = flipud(findobj(gca,'Tag','Upper Whisker'));
set(h,'linewidth',1,'linestyle','-')
for k=1:length(h)
    ydata = get(h(k),'XData');
    ydata(2) = p595z(k,2);
    set(h(k),'XData',ydata);
end
h = flipud(findobj(gca,'Tag','Upper Adjacent Value'));
set(h,'linewidth',1,'marker','none')
%h.CapSize = 12
for k=1:length(h)
    ydata = get(h(k),'XData');
    ydata(:) = p595z(k,2);
    set(h(k),'XData',ydata);
end
h = flipud(findobj(gca,'Tag','Lower Whisker'));
set(h,'linewidth',1,'linestyle','-')
for k=1:length(h)
    ydata = get(h(k),'XData');
    ydata(1) = p595z(k, 1);
    set(h(k),'XData',ydata);
end
h = flipud(findobj(gca,'Tag','Lower Adjacent Value'));
set(h,'linewidth',1,'marker','none')
for k=1:length(h)
    ydata = get(h(k),'XData');
    ydata(:) = p595z(k, 1);
    set(h(k),'XData',ydata);
end

if strcmp(opt, 'full')
    for n=1:nStrat
        % Mode add
        [~, idMainMode] = max(Pzz(idModez{n}, n)); 
        modes = locs(idModez{n});
        l2 = plot(modes, n, 'sk', 'MarkerFaceColor', [1 1 0.7], ...
                  'MarkerSize', 7, 'DisplayName','Modes');
        
        l21 = plot(modes(idMainMode), n, 'sk', 'MarkerFaceColor', 'y', ...
                   'MarkerSize', 8, 'DisplayName','Main');
        
        if mod(n, 2) == 0 % shallow fault
            krat = 3;            
        else            % medium depth fault
            krat = 10;
        end
        % Spe02 add
        SGR = mean(getSGR(thickness{n}, vcl{n}));
        logSpe02 = log10(getPermSpe02(SGR, zf(n,1), zmax{n}{1}(1)));
        l3 = plot(log10(10^(logSpe02)*krat), n, 'ok', ...
                  'MarkerFaceColor', [0.8 0.8 0.8], ...
                  'MarkerSize', 4, 'DisplayName','Spe02');
        
        % Grant19 LDL add
        A = log10(getPermSpe02(0, zf(n,1), zmax{n}{1}(1)));
        if mod(n, 2) == 0
            B = -0.15;
            c = 0;
            d = 5;
            E = -0.02;
        else
            B = -0.1;
            c = 0;
            d = 3.5;
            E = -0.012;
        end
        
        SGR = 100*SGR;
        if SGR < 15
            logGra19 = A - c + SGR*(B*15-(d-c))/15;
        elseif SGR >= 15 && SGR <= 40
            logGra19 = A - d + B*SGR;
        else
            logGra19 = A - d + 40*B + E*(SGR-40);
        end
        l4 = plot(log10(10^(logGra19)*krat), n, 'dk', 'MarkerFaceColor', [0, 138, 255]/255, ...
            'MarkerSize', 3, 'DisplayName','LDL');
    end
    %ylabel('$V_\mathrm{cl}$ [-]', 'fontSize', sz(2), latx{:})
    xlabel('$\log_{10}(k_{zz}$ [mD])', 'fontSize', sz(2), latx{:})
    %leg = legend([l1(1), l2(1), l21(1), l3, l4]);
    %set(leg.BoxFace, 'ColorType','truecoloralpha', ...
    %    'ColorData', uint8(255*[1;1;1;.5]));
    grid on
end
set(fh, 'position', [500, 500, 600, 375]);

end
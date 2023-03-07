function plotPcKr(sg, pc, krw, krn, fluids, nsim, ids, kr_mode)
%
%
%
if strcmp(kr_mode, 'CL')
    id = 2;
else
    id = 1;
end
latx = {'interpreter', 'latex'};
ct1 = 0;
cts = 0;
ctc = 0;
sw = cell(nsim, 1);
for n=1:nsim
    sw{n} = 1 - sg{n};
    if ~ids(n,1) && ~ids(n,2)
        ct1 = ct1+1;
        if ct1 == 1
            f1=figure(100);
            hold on
            p1 = plot(sg{n}, pc{n}/1e5, '--', 'color', [0.9 0.9 0.9], ...
                      'DisplayName', 'upscaled');
            f2=figure(101);
            subplot(1,2,1)
            hold on
            p21 = plot(sw{n}, krw{n}(:,id), '--', 'color', [0.9 0.9 0.9], ...
                      'DisplayName', 'upscaled');
            subplot(1,2,2)
            hold on
            p31 = plot(sw{n}, krn{n}(:,id), '--', 'color', [0.9 0.9 0.9], ...
                      'DisplayName', 'upscaled');
        else
            figure(100)
            plot(sg{n}, pc{n}/1e5, '--', 'color', [0.9 0.9 0.9])
            figure(101)
            subplot(1,2,1)
            plot(sw{n}, krw{n}(:,id), '--', 'color', [0.9 0.9 0.9])
            subplot(1,2,2)
            plot(sw{n}, krn{n}(:,id), '--', 'color', [0.9 0.9 0.9])
        end
    end
end

pcs = nan(100, nsim);
krws = nan(100, nsim);
krns = nan(100, nsim);
sgs = nan(100, nsim);
for n=1:nsim
    opt.nval = numel(pc{n});
    pcs(1:opt.nval, n) = pc{n};
    krws(1:opt.nval, n) = krw{n}(:,id);
    krns(1:opt.nval, n) = krn{n}(:,id);
    sgs(1:opt.nval, n) = sg{n};
    if ids(n,1)
        cts = cts+1;
        if cts == 1
            figure(100)
            p2 = plot(sg{n}, pc{n}/1e5, '-', 'color', [0.3 0.5 0.8], ...
                      'linewidth', 1, 'DisplayName', 'Mode high');
            figure(101)
            subplot(1,2,1)
            p22 = plot(sw{n}, krw{n}(:,id), '-', 'color', [0.3 0.5 0.8], ...
                      'linewidth', 1, 'DisplayName', 'Mode high');
            subplot(1,2,2)
            p32 = plot(sw{n}, krn{n}(:,id), '-', 'color', [0.3 0.5 0.8], ...
                      'linewidth', 1, 'DisplayName', 'Mode high');
        else
            figure(100)
            plot(sg{n}, pc{n}/1e5, '-', 'color', [0.3 0.5 0.8]);
            figure(101)
            subplot(1,2,1)
            plot(sw{n}, krw{n}(:,id), '-', 'color', [0.3 0.5 0.8]);
            subplot(1,2,2)
            plot(sw{n}, krn{n}(:,id), '-', 'color', [0.3 0.5 0.8]);
        end  
    elseif ids(n,2)
        ctc = ctc+1;
        if ctc == 1
            figure(100)
            p3 = plot(sg{n}, pc{n}/1e5, '-', 'color', [0.5 0.5 0.5], ...
                      'linewidth', 1, 'DisplayName', 'Mode low');
            figure(101)
            subplot(1,2,1)
            p23 = plot(sw{n}, krw{n}(:,id), '-', 'color', [0.5 0.5 0.5], ...
                       'linewidth', 1, 'DisplayName', 'Mode low');
            subplot(1,2,2)
            p33 = plot(sw{n}, krn{n}(:,id), '-', 'color', [0.5 0.5 0.5], ...
                       'linewidth', 1, 'DisplayName', 'Mode low');
        else
            figure(100)
            plot(sg{n}, pc{n}/1e5, '-', 'color', [0.5 0.5 0.5]);
            figure(101)
            subplot(1,2,1)
            plot(sw{n}, krw{n}(:,id), '-', 'color', [0.5 0.5 0.5]);
            subplot(1,2,2)
            plot(sw{n}, krn{n}(:,id), '-', 'color', [0.5 0.5 0.5]);
        end  
    end
end

% Pc
figure(100)
pc_avg_hi = median(pcs(:, ids(:, 1)), 2, 'omitnan');
sg_avg_hi = median(sgs(:, ids(:,1)), 2, 'omitnan');
pc_avg_lo = median(pcs(:, ids(:, 2)), 2, 'omitnan');
sg_avg_lo = median(sgs(:, ids(:, 2)), 2, 'omitnan');
%p2 = plot(sg(:,ids(1)), pc(:,ids(1))/1e5, '-', 'color', [0.5 0.5 0.5], ...
%          'linewidth', 1.5, 'DisplayName', 'Mode high');
%p3 = plot(sg(:,ids(2)), pc(:,ids(2))/1e5, '-', 'color', [0 0 0], ...
%          'linewidth', 1.5, 'DisplayName', 'Mode low');
p4 = plot(sg_avg_hi, pc_avg_hi/1e5, '-', 'color', [.1 .2 .9], ...
          'linewidth', 1.5, 'DisplayName', 'High mode');
p5 = plot(sg_avg_lo, pc_avg_lo/1e5, '-', 'color', [0 0 0], ...
          'linewidth', 1.5, 'DisplayName', 'Low mode');
colr = copper(8);
sgs = [sg{1}(1) 1e-3 0.01:0.01:sg{1}(end)];
p6 = plot(sgs, fluids{1}.pcOG{1}(sgs)/1e5, '-', ...
    'color', colr(7, :), 'linewidth', 1.5, 'DisplayName', 'Sand material example');
p7 = plot(sgs(sgs<0.58), fluids{1}.pcOG{2}(sgs(sgs<0.58))/1e5, '-', ...
    'color', colr(4, :), 'linewidth', 1.5, 'DisplayName', 'Clay smear example');
text(0.02, 10^2, ['P = ' num2str(round(ctc/nsim, 3))], latx{:}, ...
     'fontsize', 11, 'color', 'k')
text(0.55, 10^1, ['P = ' num2str(round(cts/nsim, 3))], latx{:}, ...
     'fontsize', 11, 'color', [0.3 0.5 0.8])
set(gca, 'YScale', 'log')
ylim([1e-1, 6e2])
yticks([1e-1 1e0 1e1 1e2 6e2])
xlim([0 1]) 
xticks(0:.1:1)
hold off
grid on
ylabel('$P_\mathrm{c}$ [bar]', 'fontsize', 14, latx{:})
xlabel('$S_\mathrm{g}$ [-]', 'fontsize', 14, latx{:})
h=legend([p1 p4 p5 p6 p7], 'location', 'southeast');
set(h.BoxFace, 'ColorType','truecoloralpha', ...
            'ColorData', uint8(255*[1;1;1;.5])); 
f1.Position = [100 100 400 350];

% Kr
figure(101)
krw_avg_hi = median(krws(:, ids(:, 1)), 2, 'omitnan');
krw_avg_lo = median(krws(:, ids(:, 2)), 2, 'omitnan');
krn_avg_hi = median(krns(:, ids(:, 1)), 2, 'omitnan');
krn_avg_lo = median(krns(:, ids(:, 2)), 2, 'omitnan');
sw_avg_hi =  1 - sg_avg_hi;
sw_avg_lo =  1 - sg_avg_lo;
subplot(1,2,1)
p24 = plot(sw_avg_hi, krw_avg_hi, '-', 'color', [.1 .2 .9], ...
          'linewidth', 1.5, 'DisplayName', 'High mode');
p25 = plot(sw_avg_lo, krw_avg_lo, '-', 'color', [0 0 0], ...
          'linewidth', 1.5, 'DisplayName', 'Low mode');
colr = copper(8);
sgs = [sg{1}(1) 1e-3 0.01:0.01:sg{1}(end)];
sws = 1 - sgs;
p26 = plot(sws, fluids{1}.krOG{1}(sws), '-', ...
    'color', colr(7, :), 'linewidth', 1, 'DisplayName', 'Sand material');
p27 = plot(sws(sws>0.421), fluids{1}.krOG{2}(sws(sws>0.421)), '-', ...
    'color', colr(4, :), 'linewidth', 1, 'DisplayName', 'Clay smear');
text(0.02, 0.05, ['P$_\mathrm{lo}$ = ' num2str(round(ctc/nsim, 3))], latx{:}, ...
     'fontsize', 11, 'color', 'k')
text(0.02, 0.15, ['P$_\mathrm{hi}$ = ' num2str(round(cts/nsim, 3))], latx{:}, ...
     'fontsize', 11, 'color', [0.3 0.5 0.8])
ylim([0, 1])
yticks(0:0.1:1)
xlim([0 1]) 
xticks(0:.1:1)
hold off
grid on
ylabel('$k_\mathrm{r,w}$ [-]', 'fontsize', 14, latx{:})
xlabel('$S_\mathrm{w}$ [-]', 'fontsize', 14, latx{:})
h2=legend([p21 p24 p25 p26 p27], 'location', 'southeast');
set(h2.BoxFace, 'ColorType','truecoloralpha', ...
    'ColorData', uint8(255*[1;1;1;.5])); 

subplot(1,2,2)
p34 = plot(sw_avg_hi, krn_avg_hi, '-', 'color', [.1 .2 .9], ...
          'linewidth', 1.5, 'DisplayName', 'High mode');
p35 = plot(sw_avg_lo, krn_avg_lo, '-', 'color', [0 0 0], ...
          'linewidth', 1.5, 'DisplayName', 'Low mode');
colr = copper(8);
sgs = [sg{1}(1) 1e-3 0.01:0.01:sg{1}(end)];
sws = 1 - sgs;
p36 = plot(sws, fluids{1}.krG{1}(sgs), '-', ...
    'color', colr(7, :), 'linewidth', 1.5, 'DisplayName', 'Sand material');
p37 = plot(sws(sws>0.421), fluids{1}.krG{2}(sgs(sgs<0.58)), '-', ...
    'color', colr(4, :), 'linewidth', 1.5, 'DisplayName', 'Clay smear');
ylim([0, 1])
yticks(0:0.1:1)
xlim([0 1]) 
xticks(0:.1:1)
hold off
grid on
ylabel('$k_\mathrm{r,n}$ [-]', 'fontsize', 14, latx{:})
xlabel('$S_\mathrm{w}$ [-]', 'fontsize', 14, latx{:})

f2.Position = [100 100 700 300];
end
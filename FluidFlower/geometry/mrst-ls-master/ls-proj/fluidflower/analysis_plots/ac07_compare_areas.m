%
%
%

clear, close all

% 1. Experiment
A_measured = [281.85 0 33.57 140.48 722.71 31.6 86.65 307.29];

% 2. Simulation
pth = '/home/lsalo/matlab/sim_data/mrst/fluidflower/analysis/AC07/';
mdl = {'model1/', 'model2/', 'model3/'};
fnames1 = {'areas_exp2_mesh4_modelcase1_D1e-09_I1m0.85_I2m0.85_ESF0.1mm_ESFPceSg1e-4_C_0.66mm_Cf_0.2mm_E_1.45mm_Fsupinf_3mm_Fmed_2mm_FPc0_EPc0.125_Cpc0.5_CfinfPce5mb_CfsupPce3.5mb'};
fnames2 = {'areas_exp2_mesh4_modelcase2_D1e-09_I1m0.85_I2m0.85_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc0.33_Fksupinf1.5_Fkmid1_Ek1.5_Cfk0.33_CfinfPce5mb_CfsupPce3.5mb'};
fnames3 = {'areas_exp2_mesh4_modelcase3_D1e-09_I1m0.85_I2m0.85_repPhi_repSwcSgt_model1krgkrw_Cpc1.5Sge1e-4_ESFpc2_ESFk0.33_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce5mbSge1e-3_CfsupPce3.5mbSge1e-3'};   
areas = cell(1,3);
mae = nan(1,3);
for n=1:3
    if n==1,     dirnames = fnames1;
    elseif n==2, dirnames = fnames2;  
    elseif n==3, dirnames = fnames3;   
    end
    for m=1:numel(dirnames)
        data = readmatrix([pth mdl{n} dirnames{m} '.csv']);
        areas{n}(m,:) = [data(1,2) data(2,3) data(3,3) data(4,3) ...
                         data(5,2) data(6,3) data(7,3) data(8,3)];
        mae(m,n) = sum(abs(areas{n}(m,:) - A_measured))/numel(A_measured);
    end
end

% 4. Plots
latx = {'interpreter', 'latex'};
l = numel(A_measured);
id = [1 3:l];
h=figure(17);
hold on
p1 = plot(id, A_measured(id), 'pk', 'markerSize', 12, ...
          'markerfacecolor', [0.5 0.5 0.5], 'displayname', '$E$');
p2 = plot(id, areas{1}(id), 'ob', 'markerfacecolor', 'b', ...
          'markersize', 7, 'displayname', 'm$_\mathrm{I}$');
p3 = plot(id, areas{2}(id), 'sc', 'markerfacecolor', 'c', ...
          'markersize', 6, 'displayname', 'm$_\mathrm{II}$');
p4 = plot(id, areas{3}(id), 'dm', 'markerfacecolor', 'm', ...
          'markersize', 4, 'displayname', 'm$_\mathrm{III}$');
hold off
grid on
xlim([1 l]); xticks(1:l)
ticklabels = {'$A_\mathrm{g}$ (F inf)', '$A_\mathrm{g}$ (F mid, left)', '$A_\mathrm{g}$ (F mid, right)', ...
              '$A_\mathrm{g}$ (F sup)', '$A_\mathrm{d}$ (F inf)', '$A_\mathrm{d}$ (F mid, left)', ...
              '$A_\mathrm{d}$ (F mid, right)', '$A_\mathrm{d}$ (F sup)'};
set(gca, 'XTickLabel', ticklabels, 'TickLabelInterpreter', 'latex', 'fontsize', 11);
set(gca,'YScale','log'); ylim([10 1000]); yticks([10 100 1000])
%xlabel('Region [-]', latx{:}, 'fontsize', 12)
ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
legend([p1 p2 p3 p4], 'fontsize', 10, latx{:}, 'location', 'best')
set(h, 'position', [100 100 350 325])

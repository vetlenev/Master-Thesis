%
%
%

clear, close all
inj_type = 1;

% 1. Experiment
if inj_type == 1 % AC02
    A_measured =   [91.06 3.22 108.64 211.79 33.21 251.53]; %cm2
elseif inj_type == 3 % AC07
    A_measured = [148.67 10.43 102.86 299.01 71.76 284.66];
end

% 2. Simulation
pth = 'C:/Users/lsalo/matlab/sim_data/mrst/fluidflower/results/AC02/';
mdl = {'model1/', 'model2/', 'model3/'};

fnames1 = {'areas_exp1_mesh4_modelcase1_D1e-09_I1m1_I2m1_ESF0.2mm_ESFPceSg1e-5_C_0.66mm_E_1.45mm_F_1.77mm';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m1_I2m1_ESF0.2mm_ESFPceSg1e-5_C_0.66mm_E_1.45mm_F_1.77mm_FPc0';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m1_I2m1_ESF0.2mm_ESFPceSg1e-5_C_0.66mm_E_1.45mm_F_1.77mm_FPc0.1';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-5_C_0.66mm_E_1.45mm_F_1.77mm_FPc0';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_1.77mm_FPc0_EPc0.5_CPc0.5';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_1.77mm_FPc0_EPc0.25_CPc0.5';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_1.77mm_FPc0_EPc0.125_CPc0.5';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_2.5mm_FPc0_EPc0.125_CPc0.5';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_2mm_FPc0_EPc0.125_CPc0.5';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_2.2mm_FPc0_EPc0.125_CPc0.5';
          'areas_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.1mm_ESFPceSg1e-4_C_0.66mm_Cf_0.2mm_E_1.45mm_Fsupinf_3mm_Fmed_2mm_FPc0_EPc0.125_Cpc0.5';
          'tFingers_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_1.77mm_FPc0_EPc0.125_CPc0.5';
          'tFingers_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_2mm_FPc0_EPc0.125_CPc0.5';
          'tFingers_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_2.2mm_FPc0_EPc0.125_CPc0.5';
          'tFingers_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_2.5mm_FPc0_EPc0.125_CPc0.5';
          'tFingers_exp1_mesh4_modelcase1_D1e-09_I1m0.8_I2m0.8_ESF0.1mm_ESFPceSg1e-4_C_0.66mm_Cf_0.2mm_E_1.45mm_Fsupinf_3mm_Fmed_2mm_FPc0_EPc0.125_Cpc0.5'};
      
fnames2 = {'areas_exp1_mesh4_modelcase2_D1e-09_I1m1_I2m1_repPhi_repk_repPc_model1krgkrw';
           'areas_exp1_mesh4_modelcase2_D1e-09_I1m1_I2m1_repPhi_repk_model1krgkrw_FPc0';
           'areas_exp1_mesh4_modelcase2_D1e-09_I1m0.85_I2m0.85_repPhi_repk_model1krgkrw_FPc0';
           'areas_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.5_Cpc0.5';
           'areas_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.25_Cpc0.5';
           'areas_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.125_Cpc0.33';
           'areas_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.18_Cpc0.33_Fk2';
           'areas_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc0.33_Fksupinf1.5_Fkmid1_Ek1.5_Cfk0.33';
           'tFingers_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.25_Cpc0.5';
           'tFingers_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.125_Cpc0.33';
           'tFingers_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.18_Cpc0.33_Fk2';
           'tFingers_exp1_mesh4_modelcase2_D1e-09_I1m0.9_I2m0.9_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc0.33_Fksupinf1.5_Fkmid1_Ek1.5_Cfk0.33'};
       
fnames3 = {'areas_exp1_mesh4_modelcase3_D1e-09_I1m1_I2m1_repPhi_repk_repPc_repSwcSgt_model1krgkrw';
           'areas_exp1_mesh4_modelcase3_D1e-09_I1m0.85_I2m0.85_repPhi_repk_repPc_repSwcSgt_model1krgkrw';
           'areas_exp1_mesh4_modelcase3_D1e-09_I1m0.9_I2m0.9_repPhi_repk_repPc_repSwcSgt_model1krgkrw';
           'areas_exp1_mesh4_modelcase3_D1e-09_I1m0.9_I2m0.9_repPhi_repk_repSwcSgt_model1krgkrw_Cpc1.5';
           'areas_exp1_mesh4_modelcase3_D1e-09_I1m0.9_I2m0.9_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33_Fk2';
           'areas_exp1_mesh4_modelcase3_D1e-09_I1m0.875_I2m0.875_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33_Fk1.5_Cfk0.5';
           'areas_exp1_mesh4_modelcase3_D1e-09_I1m0.875_I2m0.875_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2';
           'tFingers_exp1_mesh4_modelcase3_D1e-09_I1m1_I2m1_repPhi_repk_repPc_repSwcSgt_model1krgkrw';
           'tFingers_exp1_mesh4_modelcase3_D1e-09_I1m0.85_I2m0.85_repPhi_repk_repPc_repSwcSgt_model1krgkrw';
           'tFingers_exp1_mesh4_modelcase3_D1e-09_I1m0.9_I2m0.9_repPhi_repk_repPc_repSwcSgt_model1krgkrw';
           'tFingers_exp1_mesh4_modelcase3_D1e-09_I1m0.9_I2m0.9_repPhi_repk_repSwcSgt_model1krgkrw_Cpc1.5';
           'tFingers_exp1_mesh4_modelcase3_D1e-09_I1m0.9_I2m0.9_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33_Fk2';
           'tFingers_exp1_mesh4_modelcase3_D1e-09_I1m0.875_I2m0.875_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33_Fk1.5_Cfk0.5';
           'tFingers_exp1_mesh4_modelcase3_D1e-09_I1m0.875_I2m0.875_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2'};
tf_start = [7 5 4];
n_it = [11, 8, 7];
areas = cell(1,3);
t_fingers_s = cell(1,3);
mae = nan(11,3);
for n=1:3
    if n==1,     dirnames = fnames1;
    elseif n==2, dirnames = fnames2;  
    elseif n==3, dirnames = fnames3;   
    end
    for m=1:numel(dirnames)
        data = readmatrix([pth mdl{n} dirnames{m} '.csv']);
        if m <= n_it(n)
        	areas{n}(m,:) = [data(1,2) data(2,3) data(3,3) ...
                             data(4,2) data(5,3) data(6,3)];
            mae(m,n) = sum(abs(areas{n}(m,:) - A_measured))/numel(A_measured);
        else
            t_fingers_s{n}(m-n_it(n),:) = data(:,2);
        end
    end
end

% 4. Plots
latx = {'interpreter', 'latex'};
h=figure(17);
subplot(3,2,1)  %Ag_sup
hold on
p1 = plot([1 max(n_it)], [A_measured(3) A_measured(3)], ...
          '-k', 'linewidth', 1, 'displayname', '$E$');
p2 = plot(1:n_it(1), areas{1}(:,3), 'ob', 'markerfacecolor', ...
          'b', 'markersize', 7, 'displayname', 'm$_\mathrm{I}$');
p3 = plot(1:n_it(2), areas{2}(:,3), 'sc', 'markerfacecolor', 'c', ...
          'markersize', 6, 'displayname', 'm$_\mathrm{II}$');
p4 = plot(1:n_it(3), areas{3}(:,3), 'dm', 'markerfacecolor', 'm', 'markersize', 4, ...
          'displayname', 'm$_\mathrm{III}$');
hold off
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
ylim([0 120]); yticks([0 40 80 120])
xlabel('Iteration no. [-]', latx{:}, 'fontsize', 12)
ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
%text(3,10,'$A_\mathrm{g}$ (F sup)', 'fontsize', 12, latx{:})
legend([p1 p2 p3 p4], 'fontsize', 10, latx{:}, 'location', 'best')

subplot(3,2,2)  %Aw_sup
hold on
plot([1 max(n_it)], [A_measured(6) A_measured(6)], '-k', 'linewidth', 1);
plot(1:n_it(1), areas{1}(:,6), 'ob', 'markerfacecolor', 'b', 'markersize', 7);
plot(1:n_it(2), areas{2}(:,6), 'sc', 'markerfacecolor', 'c', 'markersize', 6);
plot(1:n_it(3), areas{3}(:,6), 'dm', 'markerfacecolor', 'm', 'markersize', 4);
hold off
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
ylim([180 340]); yticks([180 220 260 300 340])
%xlabel('Iteration no. [-]', latx{:}, 'fontsize', 12)
%ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
%text(4,318,'$A_\mathrm{d}$ (F sup)', 'fontsize', 12, latx{:})

subplot(3,2,3)  %Ag_mid
hold on
plot([1 max(n_it)], [A_measured(2) A_measured(2)], '-k', 'linewidth', 1);
plot(1:n_it(1), areas{1}(:,2), 'ob', 'markerfacecolor', 'b', 'markersize', 7);
plot(1:n_it(2), areas{2}(:,2), 'sc', 'markerfacecolor', 'c', 'markersize', 6);
plot(1:n_it(3), areas{3}(:,2), 'dm', 'markerfacecolor', 'm', 'markersize', 4);
hold off
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
%xlabel('Iteration no. [-]', latx{:}, 'fontsize', 12)
%ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
%text(4,27,'$A_\mathrm{g}$ (F mid)', 'fontsize', 12, latx{:})

subplot(3,2,4)  %Aw_mid
hold on
plot([1 max(n_it)], [A_measured(5) A_measured(5)], '-k', 'linewidth', 1);
plot(1:n_it(1), areas{1}(:,5), 'ob', 'markerfacecolor', 'b', 'markersize', 7);
plot(1:n_it(2), areas{2}(:,5), 'sc', 'markerfacecolor', 'c', 'markersize', 6);
plot(1:n_it(3), areas{3}(:,5), 'dm', 'markerfacecolor', 'm', 'markersize', 4);
hold off
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
%xlabel('Iteration no. [-]', latx{:}, 'fontsize', 12)
%ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
%text(4,72,'$A_\mathrm{d}$ (F mid)', 'fontsize', 12, latx{:})

subplot(3,2,5)  %Ag_inf
hold on
p1 = plot([1 max(n_it)], [A_measured(1) A_measured(1)], '-k', ...
          'linewidth', 1, 'displayname', '$E$');
p2 = plot(1:n_it(1), areas{1}(:,1), 'ob', 'markerfacecolor', 'b', ...
          'markersize', 7, 'displayname', 'm$_\mathrm{I}$');
p3 = plot(1:n_it(2), areas{2}(:,1), 'sc', 'markerfacecolor', 'c', ...
          'markersize', 6, 'displayname', 'm$_\mathrm{II}$');
p4 = plot(1:n_it(3), areas{3}(:,1), 'dm', 'markerfacecolor', 'm', ...
          'markersize', 4, 'displayname', 'm$_\mathrm{III}$');
hold off
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
ylim([0 120]); yticks([0 40 80 120])
%xlabel('Iteration no. [-]', latx{:}, 'fontsize', 12)
%ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
%text(3,10,'$A_\mathrm{g}$ (F inf)', 'fontsize', 12, latx{:})
%legend([p1 p2 p3 p4], 'fontsize', 10, latx{:}, 'location', 'best')

subplot(3,2,6)  %Aw_inf
hold on
plot([1 max(n_it)], [A_measured(4) A_measured(4)], '-k', 'linewidth', 1);
plot(1:n_it(1), areas{1}(:,4), 'ob', 'markerfacecolor', 'b', 'markersize', 7);
plot(1:n_it(2), areas{2}(:,4), 'sc', 'markerfacecolor', 'c', 'markersize', 6);
plot(1:n_it(3), areas{3}(:,4), 'dm', 'markerfacecolor', 'm', 'markersize', 4);
hold off
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
ylim([200 300])
%xlabel('Iteration no. [-]', latx{:}, 'fontsize', 12)
%ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
%text(4,287,'$A_\mathrm{d}$ (F inf)', 'fontsize', 12, latx{:})
set(h, 'position', [100 100 600 600])
set(h,'PaperType','A2')

h = figure(18);
hold on
p1 = plot(1:n_it(1), mae(:,1), '-ob', 'linewidth', 1, 'markerfacecolor', 'b', ...
          'markersize', 7, 'displayname', 'm$_\mathrm{I}$');
p2 = plot(1:n_it(2), mae(1:n_it(2),2), '-sc', 'linewidth', 1, 'markerfacecolor', 'c', ...
          'markersize', 6, 'displayname', 'm$_\mathrm{II}$');
p3 = plot(1:n_it(3), mae(1:n_it(3),3), '-dm', 'linewidth', 1, 'markerfacecolor', 'm', ...
          'markersize', 4, 'displayname', 'm$_\mathrm{III}$');
hold off
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
ylim([0 60]); yticks(0:10:60)
xlabel('Iteration no. [-]', latx{:}, 'fontsize', 12)
ylabel('MAE [cm$^2$]', latx{:}, 'fontsize', 12)
legend([p1 p2 p3], 'fontsize', 10, latx{:}, 'location', 'best')
set(h, 'position', [100 100 350 250])

h=figure(19);
subplot(1,2,1)
hold on
p1 = plot([1 max(n_it)], [5.25 5.25], '-k', ...
          'linewidth', 1, 'displayname', '$E$');
p2 = plot(tf_start(1):n_it(1), t_fingers_s{1}(:,1)/hour, 'ob', 'markerfacecolor', 'b', ...
          'markersize', 7, 'displayname', 'm$_\mathrm{I}$');
p3 = plot(tf_start(2):n_it(2), t_fingers_s{2}(:,1)/hour, 'sc', 'markerfacecolor', 'c', ...
          'markersize', 6, 'displayname', 'm$_\mathrm{II}$');
p4 = plot(tf_start(3):n_it(3), t_fingers_s{3}(tf_start(3):n_it(3),1)/hour, 'dm', ...
          'markerfacecolor', 'm', 'markersize', 4, 'displayname', 'm$_\mathrm{III}$');
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
xlabel('Iteration no. [-]', latx{:}, 'fontsize', 12)
ylabel('$t$ [h]', latx{:}, 'fontsize', 12)
legend([p1 p2 p3 p4], 'fontsize', 10, latx{:}, 'location', 'best')

subplot(1,2,2)
hold on
p1 = plot([1 max(n_it)], [7 7], '-k', ...
          'linewidth', 1, 'displayname', '$E$');
p2 = plot(tf_start(1):n_it(1), t_fingers_s{1}(:,2)/hour, 'ob', 'markerfacecolor', 'b', ...
          'markersize', 7, 'displayname', 'm$_\mathrm{I}$');
p3 = plot(tf_start(2):n_it(2), t_fingers_s{2}(:,2)/hour, 'sc', 'markerfacecolor', 'c', ...
          'markersize', 6, 'displayname', 'm$_\mathrm{II}$');
p4 = plot(tf_start(3):n_it(3), t_fingers_s{3}(tf_start(3):n_it(3),2)/hour, 'dm', ...
          'markerfacecolor', 'm', 'markersize', 4, 'displayname', 'm$_\mathrm{III}$');
grid on
xlim([1 max(n_it)]); xticks(1:max(n_it))
set(h, 'position', [100 100 600 250])

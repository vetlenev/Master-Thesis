% use segmentation data from experiment to compute areas with gas and brine
% with CO2 in Box A and Box B, and compare with simulation results.
clear, close all

%% 1. Experiment
pth = 'C:/Users/lsalo/matlab/sim_data/mrst/fluidflower/results/experiments/spatial_maps/';
dir_names = {'run1/';
             'run2/';
             'run4/'
             'run5/'};
fname = 'segmentation_';
timesteps = [0:300:6*3600 7*3600:3600:48*3600 54*3600:6*3600:120*3600];
ndig = 6;
dims = [150, 280];
idA1i = 90*dims(2) + 11*10+1;
idA1f = idA1i + 17*10 - 1;
idA1 = idA1i:idA1f;
id_BoxA = repmat(idA1, 60, 1) + (0:59)'*dims(2);

%id_Gexp = reshape(1:150*280, 280, 150)';
%id_BoxA = id_Gexp(91:150, 111:280);    % check we get correct cells
id_Ar = 91:150;
id_Ac = 111:280;
id_Br = 31:90;
id_Bc = 1:110;

area_A = zeros(numel(timesteps), 3); %[water, water with co2, gas] mean
area_B = zeros(numel(timesteps), 3);
std_A  = zeros(numel(timesteps), 3);
std_B  = zeros(numel(timesteps), 3);
for n=1:numel(timesteps)
    nanrun = [];
    dig = numel(num2str(timesteps(n)));
    nzer = ndig-dig;
    run1 = readmatrix([pth dir_names{1} fname ...
                       repmat('0', 1, nzer) num2str(timesteps(n)) 's.csv']);
    run4 = readmatrix([pth dir_names{3} fname ...
                       repmat('0', 1, nzer) num2str(timesteps(n)) 's.csv']);
    Ar1 = run1(id_Ar, id_Ac);
    Ar4 = run4(id_Ar, id_Ac);
    Br1 = run1(id_Br, id_Bc);
    Br4 = run4(id_Br, id_Bc);
    if n == 121     % some experimental data not present
        Ar2 = nan(numel(id_Ar), numel(id_Ac));
        Br2 = nan(numel(id_Br), numel(id_Bc));
        nanrun = 2;
    else
        run2 = readmatrix([pth dir_names{2} fname ...
                       repmat('0', 1, nzer) num2str(timesteps(n)) 's.csv']);
        Ar2 = run2(id_Ar, id_Ac);
        Br2 = run2(id_Br, id_Bc);
    end
    if ismember(n, 69:73)
        Ar5 = nan(numel(id_Ar), numel(id_Ac));
        Br5 = nan(numel(id_Br), numel(id_Bc));
        nanrun = 5;
    else
        run5 = readmatrix([pth dir_names{4} fname ...
                           repmat('0', 1, nzer) num2str(timesteps(n)) 's.csv']);
        Ar5 = run5(id_Ar, id_Ac);
        Br5 = run5(id_Br, id_Bc);
    end
    
    % Compute areas, cm2 (only need to sum since each cell is 1cm2)
    A_wat    = [sum(sum(Ar1==0)) sum(sum(Ar2==0)) sum(sum(Ar4==0)) sum(sum(Ar5==0))];
    A_watco2 = [sum(sum(Ar1==1)) sum(sum(Ar2==1)) sum(sum(Ar4==1)) sum(sum(Ar5==1))];
    A_co2    = [sum(sum(Ar1==2)) sum(sum(Ar2==2)) sum(sum(Ar4==2)) sum(sum(Ar5==2))];
    B_wat    = [sum(sum(Br1==0)) sum(sum(Br2==0)) sum(sum(Br4==0)) sum(sum(Br5==0))];
    B_watco2 = [sum(sum(Br1==1)) sum(sum(Br2==1)) sum(sum(Br4==1)) sum(sum(Br5==1))];
    B_co2    = [sum(sum(Br1==2)) sum(sum(Br2==2)) sum(sum(Br4==2)) sum(sum(Br5==2))];
    k = 1:4;
    if nanrun == 2,     k = [1 3 4];
    elseif nanrun == 5, k = [1 2 3]; 
    end
    area_A(n,:) = [mean(A_wat(k)) mean(A_watco2(k)) mean(A_co2(k))];
    area_B(n,:) = [mean(B_wat(k)) mean(B_watco2(k)) mean(B_co2(k))]; 
    if isempty(nanrun)
        assert(sum(area_A(n,:)) == numel(id_Ar)*numel(id_Ac));
        assert(sum(area_B(n,:)) == numel(id_Br)*numel(id_Bc));
    end
    std_A(n,:)  = [std(A_wat(k)) std(A_watco2(k)) std(A_co2(k))];
    std_B(n,:)  = [std(B_wat(k)) std(B_watco2(k)) std(B_co2(k))];
end

%% 2. Simulation data
pth = 'C:/Users/lsalo/matlab/sim_data/mrst/fluidflower/results/benchmark/';
dir_names = {'bmesh5_thickVar_removeS1_modelcase1_D1e-09_Im1_SgMinHyst0.02_kdF3mm_pcD0.3/';
             'bmesh5_thickVar_removeS1_modelcase2_D1e-09_Im1_SgMinHyst0.02/';
             'bmesh5_thickVar_removeS1_modelcase3_D1e-09_Im1_SgMinHyst0.02/';
             'bmesh5_thickVar_removeS1_modelcase3_D3e-09_Im1_SgMinHyst0.02/'};
fname = 'areas.csv';
m1 = readmatrix([pth dir_names{1} fname]);
m2 = readmatrix([pth dir_names{2} fname]);
m3 = readmatrix([pth dir_names{3} fname]);
m33 = readmatrix([pth dir_names{4} fname]);

%% 3. Plot
latx = {'interpreter', 'latex'};

% Box A
std_gas = [area_A(:,3) + std_A(:,3); flipud(area_A(:,3) - std_A(:,3))];
tsh = timesteps/hour;
h = figure(23);
cb = [107, 183, 250]/255;
subplot(1,3,1)  % mobile
hold on
p0 = fill([0 5 5 0 0], [0 0 1 1 0]*2e4, cb, 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', 'inj');
p1 = fill([tsh, fliplr(tsh)], std_gas, 'k', 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', '$\pm \sigma$');
p2 = plot(tsh, area_A(:,3), '-k', 'linewidth', 1.5, ...
          'displayname', '$\overline{E}$');
p3 = plot(m1(:,1)/hour, m1(:,4), '-b', 'displayname', 'm$_{\mathrm{I}}$');
p4 = plot(m2(:,1)/hour, m2(:,4), '-c', 'displayname', 'm$_{\mathrm{II}}$');
p5 = plot(m3(:,1)/hour, m3(:,4), '-m', 'displayname', 'm$_{\mathrm{III},1}$');
p6 = plot(m33(:,1)/hour, m33(:,4), '--m', 'displayname', 'm$_{\mathrm{III},3}$');
hold off
xlim([0 72]); xticks([0 5 12 24 48 72]);
ylim([0 1500]); yticks(0:500:1500)
grid on
xlabel('$t$ [h]', latx{:}, 'fontsize', 12)
ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
text(5,50,'Gaseous CO$_2$', 'fontsize', 10, latx{:})
legend([p0, p1 p2 p3 p4 p5 p6], 'fontsize', 10, latx{:})

subplot(1,3,2)  % dissolved [SUM GAS + DISSOLVED TO MAKE A CONSISTENT WITH SIM]
std_gas = [area_A(:,3) + area_A(:,2) + std_A(:,2); ...
           flipud(area_A(:,3) + area_A(:,2) - std_A(:,2))];
hold on
fill([0 5 5 0 0], [0 0 1 1 0]*2e4, cb, 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', 'inj');
fill([tsh, fliplr(tsh)], std_gas, 'k', 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', '$\pm \sigma$');
plot(tsh, area_A(:,3) + area_A(:,2), '-k', 'linewidth', 1.5, ...
          'displayname', '$\overline{E}$');
plot(m1(:,1)/hour, m1(:,3), '-b', 'displayname', 'm$_{\mathrm{I}}$');
plot(m2(:,1)/hour, m2(:,3), '-c', 'displayname', 'm$_{\mathrm{II}}$');
plot(m3(:,1)/hour, m3(:,3), '-m', 'displayname', 'm$_{\mathrm{III},1}$');
plot(m33(:,1)/hour, m33(:,3), '--m', 'displayname', 'm$_{\mathrm{III},3}$');
hold off
xlim([0 72]); xticks([0 5 12 24 48 72]);
ylim([0 8000]); yticks(0:1000:8000)
grid on
xlabel('$t$ [h]', latx{:}, 'fontsize', 12)
ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
text(5,300,'Dissolved + gaseous CO$_2$', 'fontsize', 10, latx{:})
%legend([p1 p2], 'fontsize', 10, latx{:})

subplot(1,3,3)  % water
fill([0 5 5 0 0], [0 0 1 1 0]*2e4, cb, 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', 'inj');
std_gas = [area_A(:,1) + std_A(:,1); flipud(area_A(:,1) - std_A(:,1))];
hold on
fill([tsh, fliplr(tsh)], std_gas, 'k', 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', '$\pm \sigma$');
plot(tsh, area_A(:,1), '-k', 'linewidth', 1.5, ...
          'displayname', '$\overline{E}$');
plot(m1(:,1)/hour, m1(:,2), '-b', 'displayname', 'm$_{\mathrm{I}}$');
plot(m2(:,1)/hour, m2(:,2), '-c', 'displayname', 'm$_{\mathrm{II}}$');
plot(m3(:,1)/hour, m3(:,2), '-m', 'displayname', 'm$_{\mathrm{III},1}$');
plot(m33(:,1)/hour, m33(:,2), '--m', 'displayname', 'm$_{\mathrm{III},3}$');
hold off
xlim([0 72]); xticks([0 5 12 24 48 72]);
ylim([0 11000]); yticks(0:2000:10000)
grid on
xlabel('$t$ [h]', latx{:}, 'fontsize', 12)
ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
text(5,500,'Pure water', 'fontsize', 10, latx{:})
%legend([p1 p2], 'fontsize', 10, latx{:})
set(h, 'position', [100 100 900 300])

std_gas = [area_B(:,3) + std_B(:,3); flipud(area_B(:,3) - std_B(:,3))];
h = figure(24);
subplot(1,3,1)  % mobile
hold on
p0 = fill([0 5 5 0 0], [0 0 1 1 0]*2e4, cb, 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', 'inj');
p1 = fill([tsh, fliplr(tsh)], std_gas, 'k', 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', '$\pm \sigma$');
p2 = plot(tsh, area_B(:,3), '-k', 'linewidth', 1.5, ...
          'displayname', '$\overline{E}$');
p3 = plot(m1(:,1)/hour, m1(:,7), '-b', 'displayname', 'm$_{\mathrm{I}}$');
p4 = plot(m2(:,1)/hour, m2(:,7), '-c', 'displayname', 'm$_{\mathrm{II}}$');
p5 = plot(m3(:,1)/hour, m3(:,7), '-m', 'displayname', 'm$_{\mathrm{III},1}$');
p6 = plot(m33(:,1)/hour, m33(:,7), '--m', 'displayname', 'm$_{\mathrm{III},3}$');
hold off
xlim([0 72]); xticks([0 5 12 24 48 72]);
ylim([0 70]); yticks(0:10:70)
%set(gca,'Yscale','log')
grid on
xlabel('$t$ [h]', latx{:}, 'fontsize', 12)
ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
text(12,5,'Gaseous CO$_2$', 'fontsize', 10, latx{:})
legend([p0 p1 p2 p3 p4 p5 p6], 'fontsize', 10, latx{:})

subplot(1,3,2)  % dissolved [SUM GAS + DISSOLVED TO MAKE A CONSISTENT WITH SIM]
std_gas = [area_B(:,3) + area_B(:,2) + std_B(:,2); ...
           flipud(area_B(:,3) + area_B(:,2) - std_B(:,2))];
hold on
fill([0 5 5 0 0], [0 0 1 1 0]*2e4, cb, 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', 'inj');
fill([tsh, fliplr(tsh)], std_gas, 'k', 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', '$\pm \sigma$');
plot(tsh, area_B(:,3) + area_B(:,2), '-k', 'linewidth', 1.5, ...
          'displayname', '$\overline{E}$');
plot(m1(:,1)/hour, m1(:,6), '-b', 'displayname', 'm$_{\mathrm{I}}$');
plot(m2(:,1)/hour, m2(:,6), '-c', 'displayname', 'm$_{\mathrm{II}}$');
plot(m3(:,1)/hour, m3(:,6), '-m', 'displayname', 'm$_{\mathrm{III},1}$');
plot(m33(:,1)/hour, m33(:,6), '--m', 'displayname', 'm$_{\mathrm{III},3}$');
hold off
xlim([0 72]); xticks([0 5 12 24 48 72]);
ylim([0 2000]); yticks(0:500:2000)
grid on
xlabel('$t$ [h]', latx{:}, 'fontsize', 12)
ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
text(5,200,'Dissolved + gaseous CO$_2$', 'fontsize', 10, latx{:})
%legend([p1 p2], 'fontsize', 10, latx{:})

subplot(1,3,3)  % water
fill([0 5 5 0 0], [0 0 1 1 0]*2e4, cb, 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', 'inj');
std_gas = [area_B(:,1) + std_B(:,1); flipud(area_B(:,1) - std_B(:,1))];
hold on
fill([tsh, fliplr(tsh)], std_gas, 'k', 'facealpha', 0.3, ...
          'edgecolor', 'none', 'displayname', '$\pm \sigma$');
plot(tsh, area_B(:,1), '-k', 'linewidth', 1.5, ...
          'displayname', '$\overline{E}$');
plot(m1(:,1)/hour, m1(:,5), '-b', 'displayname', 'm$_{\mathrm{I}}$');
plot(m2(:,1)/hour, m2(:,5), '-c', 'displayname', 'm$_{\mathrm{II}}$');
plot(m3(:,1)/hour, m3(:,5), '-m', 'displayname', 'm$_{\mathrm{III},1}$');
plot(m33(:,1)/hour, m33(:,5), '--m', 'displayname', 'm$_{\mathrm{III},3}$');
hold off
xlim([0 72]); xticks([0 5 12 24 48 72]);
ylim([4500 7000]); yticks(5000:500:7000)
grid on
xlabel('$t$ [h]', latx{:}, 'fontsize', 12)
ylabel('$A$ [cm$^2$]', latx{:}, 'fontsize', 12)
text(12,4700,'Pure water', 'fontsize', 10, latx{:})
%legend([p1 p2], 'fontsize', 10, latx{:})
set(h, 'position', [100 100 900 300])

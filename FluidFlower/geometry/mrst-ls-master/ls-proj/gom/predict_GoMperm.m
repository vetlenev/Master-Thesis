%% Example 0: Single stratigraphic case + analysis (3D)
% This is a complete introductory example. It shows how to load the appropriate 
% MRST modules, define the inputs according to a given faulted stratigraphy, and 
% generate the output permeability distributions. A comprehensive analysis of 
% the results is also shown. The algorithm is run in 3D mode.
% 
% We first make sure that the workspace is clean:
clear
close all force
%rng('default')          % check repeatability

%% 1. Load Required MRST Modules
% First, navigate to the mrst folder and run |startup.m|. We can then load the 
% appropriate modules for generating MRST grids and upscale the permeability:
mrstModule add mrst-gui coarsegrid upscaling incomp mpfa mimetic
mrstVerbose on     % set to on for more insight in the command window

%% 2. Define Model and Upscale Permeability

% 2.1 Mandatory input parameters
% Footwall first and hangingwall next, e.g. {[footwall, FW], [hangingwall, HW]}. 
% We need to define the layer thickness, clay content, layer dip angle, fault 
% dip angle, faulting depth, and burial depth. Further details about input parameter 
% formatting, etc can always be checked from the documentation in the classes 
% and functions.   
thickness = {[35.8537 35.8537 35.8537 35.8537], [37.4901 35.2847 35.3553 35.2847]};  % [m]
vcl       = {[0.2878, 0.7, 0.2878, 0.7], [0.2878, 0.7, 0.2878, 0.7]};                 % fraction [-]
dip       = [0, -5.2221];                                           % [deg.]
faultDip  = 45.0685;                                           % [deg.]
zf        = [200, 200];                                        % [FW, HW], [m]
zmax      = {[1538.8 1513.82 1488.75 1463.99], [1538.8 1513.82 1488.75 1463.99]};    % {FW, HW}
dim       = 3;                    % dimensions (2 = 2D, 3 = 3D)
unit_plot = 'm';

% 2.2 Optional input parameters
% In this case, we indicate a maximum fault material permeability of and a correlation 
% coefficient for dependent variables:
maxPerm = [];                   % [mD]
rho     = 0.6;                  % Corr. coeff. for multivariate distributions

% 2.3 Flow upscaling options and number of simulations
U.method          = 'tpfa';     % 'tpfa'for 3D
U.coarseDims      = [1 1 1];    % desired n cells [x, y, z] in coarse grid
U.flexible        = true;       % default true, much faster but U.coarseDims
                                % will be modified in some realizations
                                % unless U.coarseDims = [1 1 1] (do not set
                                % to false in that case).
Nsim              = 1000;       % Number of 3D simulations/realizations

% 2.4 Define Stratigraphy and FaultedSection objects
% Organize the input parameters in HW and FW, and use that info to create a 
% FaultedSection object which contains all required information.
% FW and HW
footwall = Stratigraphy(thickness{1}, vcl{1}, 'Dip', dip(1), ...
                        'DepthFaulting', zf(1), 'DepthBurial', zmax{1});
hangingwall = Stratigraphy(thickness{2}, vcl{2}, 'Dip', dip(2), 'IsHW', 1, ...
                           'NumLayersFW', footwall.NumLayers, ...
                           'DepthFaulting', zf(2), 'DepthBurial', zmax{2});

% Instantiate FaultedSection object (Strati in Faulted Section)
mySect = FaultedSection(footwall, hangingwall, faultDip, 'maxPerm', maxPerm);

% 2.5 Get material distributions
% We use the inputs to constrain the ranges and distributions for each of the 
% intermediate variables.
% Get material property distributions
mySect = mySect.getMatPropDistr();
% Get along-strike segmentation
nSeg = getNSeg(mySect.Vcl, mySect.IsClayVcl, mySect.DepthFaulting);

% 2.6 Generate intermediate variable samples, calculate smear dimensions 
%     and upscale permeability.
% We create two container variables (faults and smears) where we'll save all 
% data for each realization. For each realization, the code defines a Fault object, 
% generates intermediate variable samples, calculates the smear dimensions, and, 
% within upscaleSmearPerm, generates a fault material distribution consistent 
% with the inputs and upscales the permeability.
% Generate fault object with properties for each realization
assert(dim==3);
faultSections = cell(Nsim, 1);
smears = cell(Nsim, 1);
faults = cell(Nsim, 1);
Us = cell(Nsim, 1);
nSeg_fcn = nSeg.fcn;
U_flex = U.flexible;
tstart = tic;
parfor n=1:Nsim    % parfor allowed if you have the parallel computing toolbox
%for n=1
    % Instantiate fault section and get segmentation for this realization
    myFaultSection = Fault2D(mySect, faultDip);
    myFault = Fault3D(myFaultSection, mySect);
    if U_flex
        [myFault, Us{n}] = myFault.getSegmentationLength(U, nSeg_fcn);
    else
        myFault = myFault.getSegmentationLength(U, nSeg_fcn);
    end
    G = [];
    for k=1:numel(myFault.SegLen)
        % Get material property (intermediate variable) samples, and fix
        % fault thickness of current realization (3D only).
        myFaultSection = myFaultSection.getMaterialProperties(mySect, 'corrCoef', rho);
        myFaultSection.MatProps.thick = myFault.Thick;
        if isempty(G)
            G = makeFaultGrid(myFault.Thick, myFault.Disp, ...
                              myFault.Length, myFault.SegLen, U);
        end
        
        % Generate smear object with T, Tap, L, Lmax
        smear = Smear(mySect, myFaultSection, G, 1);
        
        % Place fault materials and assign cell-based properties in 2D section
        myFaultSection = myFaultSection.placeMaterials(mySect, smear, G);
        
        % Extrude 2D section to fill current segment
        myFault = myFault.assignExtrudedVals(G, myFaultSection, k);
        
        % Save results
        faultSections{n}{k} = myFaultSection;
        smears{n}{k} = smear;
    end

    % Compute 3D upscaled permeability distribution
    if U_flex
        myFault = myFault.upscaleProps(G, Us{n});
    else
        myFault = myFault.upscaleProps(G, U);
    end
    
    % Save results
    faults{n} = myFault;
    if mod(n, 50) == 0
        disp(['Simulation ' num2str(n) ' / ' num2str(Nsim) ' completed.'])
    end
end
telapsed = toc(tstart);

%% 3. Output Analysis
% 3.1 Visualize stratigraphy and fault (with thickness corresponding to 1st realization)
mySect.plotStrati(faults{1}.Thick, faultDip, unit_plot);  

% 3.2 Visualize intermediate variables
% We define a given parent material (id from 1 to n of materials in stratigraphy), 
% and generate histograms and correlation matrix plots.
layerId = 4;                                            
plotMatPropsHist(faultSections, smears, mySect, layerId, dim) 
% MatProps correlations
[R, P] = plotMatPropsCorr(faultSections, mySect, layerId, dim);
if dim==3
    plotSeg(faults, nSeg)
end

% 3.3 Visualize fault materials
% Visualization for one realization. Choice can be 'randm' (random), 'maxX' 
% (realization with maximum upscaled permeability in across the fault), 'minX', 
% 'maxZ' or 'minZ'.
% General fault materials and perm view
plotId = selectSimId('randm', faults, Nsim);                % simulation index
%plotId = 1;
if U_flex
    faults{plotId}.plotMaterials(faultSections{plotId}{1}, mySect, ...
                                 unit_plot, Us{plotId}) 
else
    faults{plotId}.plotMaterials(faultSections{plotId}{1}, mySect, unit_plot, U) 
end

% 3.4. Visualize upscaled permeability
% Plot upscaled permeability distributions (all simulations)
plotUpscaledPerm(faults, dim, 'histonly')

% Make video of fault core materials
close all
latx = {'interpreter','latex'};
dir = 'C:/Users/lsalo/matlab/vids/gomFaultLeakage/famp5_faultMaterialExamples/';
idf = randi(Nsim, 1, 30);
v = VideoWriter([dir 'faultMaterialExamples.mp4'], 'MPEG-4');
v.Quality = 100;
v.FrameRate = 1;
open(v);
disp('Generating video...')
for k=1:numel(idf)
   % Get Grid
    obj = faults{idf(k)};
    G = makeFaultGrid(obj.Thick, obj.Disp, obj.Length, obj.SegLen, U);
    % Make plot
    hf = figure(27);
    tiledlayout(2, 1, 'Padding', 'tight', 'TileSpacing', 'tight');
    % 1. clay smears
    nexttile
    cmap = copper;
    set(gca, 'colormap', cmap(1:128, :))
    plotToolbar(G, log10(obj.Grid.perm(:,1)/(milli*darcy)), ...
                obj.Grid.isSmear, 'EdgeColor', [0.8 0.8 0.8], ...
                'EdgeAlpha', 0.1);
    xlim([0 obj.Thick]); zlim([0 obj.Disp]); ylim([0 obj.Length]);
    c = colorbar;
    %caxis([min(log10(obj.Grid.perm(:,1)/(milli*darcy))) ...
    %       max(log10(obj.Grid.perm(:,4)/(milli*darcy)))]);
    c.Label.Interpreter = 'latex';
    c.Label.String = '$\log_{10} k_{xx}$ [mD] (Clay smears)';
    c.Label.FontSize = 12;
    set(gca,'fontSize', 10)
    val = zeros(1,3);
    units = cell(1,3);
    form = cell(1,3);
    for n=1:3
        val(n) = obj.Perm(n)/(milli*darcy);
        if val(n) < 1e-3
            val(n) = val(n)*1000;
            units{n} = ' [$\mu$D]';
            form{n} = ' %1.3f';
        elseif val(n) > 999
            val(n) = val(n)/1000;
            units{n} = ' [D]';
            form{n} = ' %1.2f';
        else
            units{n} = ' [mD]';
            form{n} = ' %3.3f';
        end
    end
    xlabel('$x$ [m]', latx{:}); 
    ylabel('$y$ [m]', latx{:})
    zlabel('$z$ [m]', latx{:})
    ax = gca;
    ax.DataAspectRatio = [0.15 1 1];
    ax.ZDir = 'normal';
    view([30 20])
    title(['$k_{jj} =$ ' num2str(val(1), form{1}) units{1}, ...
           ' $\vert$ ' num2str(val(2), form{2}) units{2}, ...
           ' $\vert$ ' num2str(val(3), form{3}) units{3}], latx{:}, ...
        'fontSize', 12);
    grid on
    zlim([0 150])
    % 2. sand smears
    nexttile
    cmap = copper;
    set(gca, 'colormap', cmap(156:end, :))
    plotToolbar(G, log10(obj.Grid.perm(:,1)/(milli*darcy)), ...
                ~obj.Grid.isSmear, 'EdgeColor', [0.2 0.2 0.2], ...
        'EdgeAlpha', 0.1);
    xlim([0 obj.Thick]); zlim([0 obj.Disp]); ylim([0 obj.Length]);
    c = colorbar;
    %caxis([min(log10(obj.Grid.perm(:,1)/(milli*darcy))) ...
    %       max(log10(obj.Grid.perm(:,4)/(milli*darcy)))]);
    c.Label.Interpreter = 'latex';
    c.Label.String = '$\log_{10} k_{xx}$ [mD] (Sand smears)';
    c.Label.FontSize = 12;
    set(gca,'fontSize', 10)
    ax = gca;
    ax.DataAspectRatio = [0.15 1 1];
    ax.ZDir = 'normal';
    view([30 20])
    grid on
    zlim([0 150])
    xticks([]), yticks([]), zticks([])
    set(hf, 'units', 'pixels', 'position', [200, 20, 400, 800]);
    set(hf,'color','w'); 
    
    % Export graphics
    %exportgraphics(hf,[dir num2str(k) '.jpg'],'ContentType','image',...
    %               'Resolution', 300, 'BackgroundColor','w')
    print(hf, [dir num2str(k)], '-djpeg', '-r300')
               
    % Read graphics and save frame
    close(hf)
    I = imread([dir num2str(k) '.jpg']);
    
    % Add frame to video
    %frame = getframe(gcf);
    frame = im2frame(I);
    writeVideo(v,frame);
    close(gcf)
    
    if mod(k, 5) == 0
        disp([num2str(k) '/' num2str(numel(idf)) ' frames completed.'])
    end
end
close(v)
disp(['Done. Video file saved to ' dir]);   

% Fit probability distribution to permeability 
perms = cell2mat(cellfun(@(x) x.Perm, faults, ...
                 'UniformOutput', false)) ./ (milli*darcy);
perm_log = log10(perms);
%pd = fitdist(perm_log(:,1), 'Kernel', 'kernel', 'normal', 'width', 0.1); 
pdx = fitdist(perm_log(:,1), 'Kernel', 'kernel', 'normal', 'width', 0.08); 
pdy = fitdist(perm_log(:,2), 'Kernel', 'kernel', 'normal'); 
pdz = fitdist(perm_log(:,3), 'Kernel', 'kernel', 'normal', 'width', 0.1); 
k_rnd = [random(pdx, [1000,1]), random(pdy, [1000,1]), ...
         random(pdz, [1000,1])];

% Compare figure
% Histogram parameters
nbins = 50;
logMinP = -6; %min(min(K));
logMaxP = 2; %max(max(K));
edges = linspace(fix(logMinP)-1, fix(logMaxP)+1, nbins);
labls = ["$\log_{10}(k_{xx}$ [mD])", ...
        "$\log_{10}(k_{yy}$ [mD])", ...
        "$\log_{10}(k_{zz}$ [mD])"];
latx = {'interpreter','latex'};

% prepare prob dist
x = logMinP:.1:logMaxP;
pdx_plot = pdf(pdx, x)*diff(edges(1:2));         % P = pdf*binwidth
pdy_plot = pdf(pdy, x)*diff(edges(1:2));
pdz_plot = pdf(pdz, x)*diff(edges(1:2));

figure(1238)
% x
subplot(3,2,1)
p1 = histogram(perm_log(:,1), edges, 'Normalization', 'probability', ...
          'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1, 'displayname', 'PREDICT');
hold on
p2 = plot(x,pdx_plot,'-g', 'linewidth', 1, 'displayname', 'fit');
hold off
xlabel(labls(1), latx{:}, 'fontSize', 12)
ylabel('P [-]', latx{:}, 'fontSize', 12)
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:1)
grid on
legend([p1 p2], 'fontsize', 10)
subplot(3,2,2)
hold on
p1 = histogram(perm_log(:,1), edges, 'Normalization', 'probability', ...
          'FaceColor', 'none', 'edgecolor', 'k', 'linewidth', 1.25, ...
          'displayname', 'PREDICT');
p2 = histogram(k_rnd(:,1), edges, 'Normalization', 'probability', ...
          'FaceColor', 'none', 'edgecolor', 'g', 'displayname', 'From fit');
hold off
xlabel(labls(1), latx{:}, 'fontSize', 12)
ylabel('P [-]', latx{:}, 'fontSize', 12)
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:1)
grid on
legend([p1 p2], 'fontsize', 10)
% y
rr = [255, 125, 125]/255;
subplot(3,2,3)
p1 = histogram(perm_log(:,2), edges, 'Normalization', 'probability', ...
          'FaceColor', rr, 'FaceAlpha', 1, 'displayname', 'PREDICT');
hold on
p2 = plot(x,pdy_plot,'-g', 'linewidth', 1, 'displayname', 'fit');
hold off
xlabel(labls(2), latx{:}, 'fontSize', 12)
ylabel('P [-]', latx{:}, 'fontSize', 12)
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:1)
grid on
%legend([p1 p2], 'fontsize', 10)
subplot(3,2,4)
hold on
p1 = histogram(perm_log(:,2), edges, 'Normalization', 'probability', ...
          'FaceColor', 'none', 'edgecolor', 'r', 'linewidth', 1.25, ...
          'displayname', 'PREDICT');
p2 = histogram(k_rnd(:,2), edges, 'Normalization', 'probability', ...
          'FaceColor', 'none', 'edgecolor', 'g', 'displayname', 'From fit');
hold off
xlabel(labls(2), latx{:}, 'fontSize', 12)
ylabel('P [-]', latx{:}, 'fontSize', 12)
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:1)
grid on
%legend([p1 p2], 'fontsize', 10)
% z
bb = [125, 125, 255]/255;
subplot(3,2,5)
p1 = histogram(perm_log(:,3), edges, 'Normalization', 'probability', ...
          'FaceColor', bb, 'FaceAlpha', 1, 'displayname', 'PREDICT');
hold on
p2 = plot(x,pdz_plot,'-g', 'linewidth', 1, 'displayname', 'fit');
hold off
xlabel(labls(3), latx{:}, 'fontSize', 12)
ylabel('P [-]', latx{:}, 'fontSize', 12)
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:1)
grid on
%legend([p1 p2], 'fontsize', 10)
subplot(3,2,6)
hold on
p1 = histogram(perm_log(:,3), edges, 'Normalization', 'probability', ...
          'FaceColor', 'none', 'edgecolor', 'b', 'linewidth', 1.25, ...
          'displayname', 'PREDICT');
p2 = histogram(k_rnd(:,3), edges, 'Normalization', 'probability', ...
          'FaceColor', 'none', 'edgecolor', 'g', 'displayname', 'From fit');
hold off
xlabel(labls(3), latx{:}, 'fontSize', 12)
ylabel('P [-]', latx{:}, 'fontSize', 12)
xlim([fix(logMinP)-1 fix(logMaxP)+1])
ylim([0 0.5]); yticks(0:.1:1)
grid on
%legend([p1 p2], 'fontsize', 10)

% save data
%fnme = 'famp6_kfit';
%save([fnme '.mat'], 'pdx', 'pdy', 'pdz');


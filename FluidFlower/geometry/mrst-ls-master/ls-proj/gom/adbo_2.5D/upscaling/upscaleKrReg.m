function [krw_out, krn_out, sg_ups, vpar_fit, mae_min] = upscaleKrReg(G, CG, ...
                    rock, fluid, fault3D, sg, pc, sgmax, faultSection, ...
                                                   sect, U, opt, makeplots)
%
% Procedure from appendix B, Li & Benson AWR (2015)
%

if strcmp(opt.fault, 'test')
    if strcmp(opt.dir, 'z'), id_dim = 2; idkG = 1; end
else
    if strcmp(opt.dir, 'x'),        id_dim = 1;  idkG = 1;
    elseif strcmp(opt.dir, 'y'),    id_dim = 2;  idkG = 4;
    elseif strcmp(opt.dir, 'z'),    id_dim = 3;  idkG = 6;
    end
end
isSmear = fault3D.Grid.isSmear;

%---------------------------------- VL -----------------------------------
if strcmp(opt.kr_mode, 'VL') || strcmp(opt.kr_mode, 'both')    % viscous limit
    % Get pressures for later computation of viscosities. Pressures are
    % from GoM model, we get the initial pressures for cells with centroids
    % around (+-1m) z_p depth, and average them.
    %z_p = mean([opt.zmax{1}, opt.zmax{2}]);
    if strcmp(opt.window, 'famp1')
        p = 19934555.3835938;  % [Pa]
    end
    
    % kr Regions and fractional flow inverse function
    % We consider 2 regions for kr: sand and clay (without regard to the
    % variations in Vcl of each one)
    fluid.fwInv = cell(1, 2);
    nS = pow2(5);
    sgmax_fw = zeros(1,2);
    for n=1:2
        if strcmp(opt.sg, 'sandClay')
            if n == 1   % sands
                sgmin = fluid.krPts.g(1,1);
                sgmax_fw(n) = fluid.krPts.g(1,3);
            else
                sgmin = fluid.krPts.g(2,1);
                sgmax_fw(n) = fluid.krPts.g(2,3);
            end
        end
        sgvals = linspace(sgmin, sgmax_fw(n), nS-1)';
        sgvals = [sgvals(1); sgvals(1)+1e-3; sgvals(2:end)];
        swv = 1 - sgvals;
        krnv = fluid.krG{n}(sgvals);    % kr, nonwetting phase
        krwv = fluid.krOG{n}(swv);      % wetting phase
        lnv = krnv / fluid.muG(p);      % mobility
        lwv = krwv / fluid.muO(p);
        fwv = lwv ./ (lwv + lnv);       % fractional flow
        fluid.fwInv{n} = @(fw) interp1(fwv, swv, fw);   % get Sw from fw
        fluid.fw{n} = @(Sw) interp1(swv, fwv, Sw);      % get fw from Sw
    end
    
    % Upscale relative permeability in the vertical (z) direction
    fw = [1e-7 1e-6 1e-5 1e-4 1e-3:1e-3:1e-2 2e-2:.01:.99 .999 .99999 1]';
    swc = nan(G.cells.num, 1);
    mun = fluid.muG(p);
    muw = fluid.muO(p);
    volume=sum(G.cells.volumes.*rock.poro);
    sw = zeros(numel(fw), 1);
    krw = zeros(numel(fw), 1);
    krn = zeros(numel(fw), 1);
    lt = zeros(numel(fw), 1);
    for n=1:numel(fw)
        % Sw, Kr and mobilities for each cell.
        swv = [fluid.fwInv{1}(fw(n)) fluid.fwInv{2}(fw(n))]; % [sand, clay]
        sgv = 1-swv;
        swc(~isSmear) = swv(1);
        swc(isSmear) = swv(2);
        sw(n) = sum(swc.*G.cells.volumes.*rock.poro)/volume;    % PV weighted avg Sw 
        krwv = [fluid.krOG{1}(swv(1)) fluid.krOG{2}(swv(2))];
        krnv = [fluid.krG{1}(sgv(1)) fluid.krG{2}(sgv(2))];
        lnv = krnv ./ mun;
        lwv = krwv ./ muw;
        ltv = lnv + lwv;   % total mobility (without k)
        
        % Single-phase steady state simulation with viscosity = 1 and
        % permeability from total mobility previously calculated (lwt).
        % Vertical (z) direction only
        rock_ups.perm(~isSmear, 1) = ltv(1);
        rock_ups.perm(isSmear, 1) = ltv(2);
        rock_ups.perm = rock_ups.perm.*fault3D.Grid.perm;    % cell total mobility (with k)
        %assert(all(rock_ups.perm > 0))
        lt(n) = myUpscalePermDim(G, CG, rock_ups, 'dim', id_dim); 
        % b = myUpscalePerm(G, CG, rock_ups, 'method', 'tpfa'); just to check
        clear rock_ups

        % Average water saturation and upscaled relative permeabilities
        krw(n) = muw*lt(n)*fw(n)/fault3D.Perm(id_dim);
        krn(n) = mun*lt(n)*(1-fw(n))/fault3D.Perm(id_dim);  
%          krw(n) = muw*lt(n)*fw(n);
%          krn(n) = mun*lt(n)*(1-fw(n));
        if krw(n) > 1
            krw(n) = 1;
        end
        %plot(sw(1:n),krw(1:n)); hold on; plot(sw(1:n),krn(1:n));
        %grid on
    end 
    assert(max(krw) <= 1)
    assert(min(krw) >= 0)
    assert(max(krn) <= 1)
    assert(min(krn) >= 0)
    assert(max(sw) <= 1)
    assert(min(sw) >= min(1-sgmax_fw)) 
    sg_ups = flipud(1 - sw);
    krw_out(:,1) = interp1(sg_ups,flipud(krw),sg);
    krn_out(:,1) = interp1(sg_ups,flipud(krn),sg);
end   

%---------------------------------- CL -----------------------------------
if strcmp(opt.kr_mode, 'CL') || strcmp(opt.kr_mode, 'both')
    % N of regions and cells in region
    reg = rock.regions.saturation;
    id_reg = unique(rock.regions.saturation);
    nreg = numel(id_reg);
    ncreg = sum(reg(:)==(1:max(reg)));
    ncreg(ncreg==0) = [];
    % Loop params
    nS = pow2(7);
    pc_val = logspace(log10(pc(2)), log10(pc(end)), nS-1);
    pc_val = [0 0.98*pc_val(1) pc_val];
    volume=sum(G.cells.volumes.*rock.poro);
    rock_ups.perm = zeros(G.cells.num, 1);
    krw = nan(nS, 1);
    krn = nan(nS, 1);
    sw = nan(nS, 1);
    for i=1:numel(pc_val)
        sg_cells = nan(G.cells.num, 1);
        krwv = nan(G.cells.num, 1);
        krnv = nan(G.cells.num, 1);
        for n=1:nreg    % different Pc region for each unit in fault
            idcell = reg == id_reg(n);
            if pc_val(i) >= fluid.pcOG{id_reg(n)}(sgmax(n))     % pcInv returns NaN
                %pcv = fluid.pcOG{id_reg(n)}(sgmax(n));
                sg_cells(idcell) = 0.999*sgmax(n);
            else
                pcv = pc_val(i);
                sg_cells(idcell) = fluid.pcInv{id_reg(n)}(pcv*ones(ncreg(n), 1));
                assert(all(sg_cells(idcell)<sgmax(n)))
            end
        end
        % Set min sg > 0 to prevent kn=NaN during single phase sims below
        if i > 1
            if any(sg_cells < 1e-5)
                sg_cells(sg_cells < 1e-5) = 1e-5;
            end
        end

        % Sw distribution at cell level
        sw_cells = 1 - sg_cells;
        
        % Pore-volume weighted Sw corresponding to pc_val(i)
        sw(i)=sum(sw_cells.*G.cells.volumes.*rock.poro)/volume;
        
        % Cell relative and phase permeabilities. Only 2 kr regions (sand
        % and clay smear)
        krwv(~isSmear) = fluid.krOG{1}(sw_cells(~isSmear));
        krwv(isSmear) = fluid.krOG{2}(sw_cells(isSmear));
        %if i > 0, krwv(krwv<1e-3) = 1e-3; end
        krnv(~isSmear) = fluid.krG{1}(sg_cells(~isSmear));
        krnv(isSmear) = fluid.krG{2}(sg_cells(isSmear)); 
        if i > 0, krnv(krnv<1e-3) = 1e-3; end
        kwv = krwv.*fault3D.Grid.perm;
        knv = krnv.*fault3D.Grid.perm;
        %kwv = krwv;
        %knv = krnv;
        
        % Calculate effective (upscaled) phase permeabilities
        rock_ups.perm = kwv;
%        assert(all(rock_ups.perm > 0))
        kw = myUpscalePermDim(G, CG, rock_ups, 'dim', id_dim);
        clear rock_ups 
        if i == 1, kn = 0;
        else
            rock_ups.perm = knv;
%            assert(all(rock_ups.perm >= 0))
            kn = myUpscalePermDim(G, CG, rock_ups, 'dim', id_dim);
            clear rock_ups
        end
        
        % Upscaled relative perm
        krw(i) = kw / fault3D.Perm(id_dim);
        % = kw;
        if krw(i) > 1
            krw(i) = 1;
        end
        krn(i) = kn / fault3D.Perm(id_dim);
        %krn(i) = kn;
    end
    assert(~any(isnan(sw)))
    assert(~any(isnan(krw)))
    assert(~any(isnan(krn)))
    assert(max(krw) <= 1)
    assert(min(krw) >= 0)
    assert(max(krn) <= 1)
    assert(min(krn) >= 0)
    assert(max(sw) <= 1)
    assert(min(sw) >= min(1-sgmax)) 
    sg_ups = 1 - sw;
    krw_out(:,2) = interp1(sg_ups,krw,sg);
    krn_out(:,2) = interp1(sg_ups,krn,sg);
end

%-------------------------------- dynamic ---------------------------------
if strcmp(opt.kr_mode, 'dyn')
    % Run dynamic simulation at high flow rate
    
    % 1. Assumptions / options
    gravity reset off
    opt.dyn_ispc = 0;           % 0: ignore pc (VL conditions); 1: consider pc
    opt.dyn_incomp_run = true;  % true: set very low compressibility for fluids and rock
    opt.dyn_perm_case = 'sand'; % none (sand + clay), geomean or sand
    opt.dyn_mrate = 1;
    opt.dyn_tsim_year = 1;
    states_plots = 0;   % plot states
    
    % Run 3D simulation
    ts = [1*day (2:2:30)*day (60:30:360)*day 1*year];
    [states, G2, rock2, fluid, state0, rate] = dynamic3Drun(G, rock , fluid, ...
                                                            ts, fault3D, opt, states_plots);
                                                        
    % Get saturations
    ncx = prod(G2.cartDims(1:2));
    ncz = G2.cartDims(3);
    %id1 = ones(ncx, 1)*(1:G2.cartDims(3));
    %idG2to1D = flipud(reshape(id1, G2.cells.num, 1));
    sz = zeros(ncz, numel(states));
    for n=1:numel(states)
        sz(:,n) = mean(reshape(states{n}.s(:,2), ncx, ncz));
    end
    zcor = linspace(min(G2.cells.centroids(:,3)), ...
        max(G2.cells.centroids(:,3)), G2.cartDims(end))';
    
    % Plot fault materials for reference
    %fault3D.plotMaterials(faultSection{1}, sect, 'm', U) 
    
    % 1D BL simulations (see buckleyLeverett1D.m example)
    [s1D, vpar] = dynamicBLrun(G2, rock2, fluid, state0, rate, ts, opt, states_plots);
    mae_mean = nan(numel(s1D), 1);
    for n=1:numel(s1D)
       mae_mean(n) = mean(sum(abs(s1D{n}-sz))./size(sz, 1)); % compute mae at each timestep and then mean
    end
    [mae_min, id_min] = min(mae_mean);
    vpar_fit = vpar(id_min, :);
    [~, id_sort] = sort(mae_mean);
    vpar_sort = vpar(id_sort(1:10), :);
    
    % Plot saturation profiles at select times
    figure(randi(1000,1))
    latx = {'interpreter', 'latex'};
    hold on
    clrs = turbo(8);
    plot(zcor, sz(:,1), '.', 'color', clrs(1,:))
    plot(zcor, sz(:,5), '.', 'color', clrs(2,:))
    plot(zcor, sz(:,10), '.', 'color', clrs(3,:))
    plot(zcor, sz(:,16), '.', 'color', clrs(4,:))
    plot(zcor, sz(:,17), '.', 'color', clrs(5,:))
    plot(zcor, sz(:,19), '.', 'color', clrs(6,:))
    plot(zcor, sz(:,21), '.', 'color', clrs(7,:))
    plot(zcor, sz(:,28), '.', 'color', clrs(8,:))

    plot(zcor, s1D{id_min}(:,1), '-', 'color', clrs(1,:))
    plot(zcor, s1D{id_min}(:,5), '-', 'color', clrs(2,:))
    plot(zcor, s1D{id_min}(:,10), '-', 'color', clrs(3,:))
    plot(zcor, s1D{id_min}(:,16), '-', 'color', clrs(4,:))
    plot(zcor, s1D{id_min}(:,17), '-', 'color', clrs(5,:))
    plot(zcor, s1D{id_min}(:,19), '-', 'color', clrs(6,:))
    plot(zcor, s1D{id_min}(:,21), '-', 'color', clrs(7,:))
    plot(zcor, s1D{id_min}(:,28), '-', 'color', clrs(8,:))
    hold off
    grid on
    xlabel('$z$ [m]', 'fontsize', 14, latx{:})
    ylabel('$S_\mathrm{g}$ [-]', 'fontsize', 14, latx{:})
    xlim([min(zcor) max(zcor)])
    ylim([0 1])
    set(gca,'xdir','reverse')
    
    % saturations and kr out
    sg_ups = linspace(0,1-vpar_fit(3),2^7);
    %sg_ups = sg;
    sw_ups = 1-sg_ups;
    sw_ups_n = (sw_ups - min(sw_ups))/(1 - min(sw_ups));
    sg_ups_n = 1 - sw_ups_n;
    krw_out = sw_ups_n.^vpar_fit(1);
    krn_out = sg_ups_n.^vpar_fit(2);
end

if makeplots == 1
    % Relperm comparison
    latx = {'interpreter', 'latex'};
    sw_ups = 1 - sg;
    h = figure(11);
    subplot(1,2,1)
    hold on
    colrs = copper(8);
    sgp = 0:.01:sgmax(1); sgp2 = 0:.01:sgmax(2);
    swp = 1 - sgp; swp2 = 1 - sgp2;
    plot(swp,fluid.krG{1}(sgp), '-', 'color', colrs(8,:), ...
         'DisplayName', '$k_{r,g}$ (sand)', 'linewidth', 1)
    plot(swp2,fluid.krG{2}(sgp2), '-', 'color', colrs(3,:), ...
         'DisplayName', '$k_{r,g}$ (clay)', 'linewidth', 1)
    if strcmp(opt.kr_mode, 'VL') || strcmp(opt.kr_mode, 'both')
        plot(sw_ups, krn_out(:,1), '-r', 'linewidth', 1.5, ...
             'DisplayName', '$k_{r,n}$ (VL)')
    end
    if strcmp(opt.kr_mode, 'CL') || strcmp(opt.kr_mode, 'both')
        plot(sw_ups, krn_out(:,2), '--', 'color', [0.5 0.1 0.1], ...
             'linewidth', 1.5, 'DisplayName', '$k_{r,n}$ (CL)')
    end
    hold off
    grid on
    ylabel('$k_{r,g}$', latx{:})
    xticks(0:.2:1)
    yticks(0:.2:1)
    xlabel('$S_\mathrm{w}$', latx{:})
    ylim([0 1])
    xlim([0 1])
    legend(latx{:})
    subplot(1,2,2)
    hold on
    plot(swp,fluid.krOG{1}(swp), '-', 'color', colrs(8,:), ...
         'DisplayName', '$k_{r,w}$ (sand)', 'linewidth', 1)
    plot(swp2,fluid.krOG{2}(swp2), '-', 'color', colrs(3,:), ...
         'DisplayName', '$k_{r,w}$ (clay)', 'linewidth', 1)
    if strcmp(opt.kr_mode, 'VL') || strcmp(opt.kr_mode, 'both')
        plot(sw_ups, krw_out(:,1), '-b', 'linewidth', 1.5, ...
             'DisplayName', '$k_{r,w}$ (VL)')
    end
    if strcmp(opt.kr_mode, 'CL') || strcmp(opt.kr_mode, 'both')
        plot(sw_ups, krw_out(:,2), '--c', 'linewidth', 1.5, ...
             'DisplayName', '$k_{r,w}$ (CL)')
    end
    hold off
    grid on
    ylabel('$k_{r,w}$', latx{:})
    xlabel('$S_\mathrm{w}$', latx{:})
    xticks(0:.2:1)
    yticks(0:.2:1)
    xlim([0 1])
    legend(latx{:})  
end

end
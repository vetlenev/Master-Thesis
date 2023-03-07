function [plts, opt, m] = bSimOpts(inj_case, mesh_size, minj, removeS, ...
                                   D, sgmin, dF, pcD, rampdown, flatMesh, ...
                                   temp_C, modelParams, batch)
%
%
%

%% Plots
plts.fig3D = @() figure('Position', [0, 0, 1300, 650]);
plts.setAxProps = @(ax) set(ax, 'View'              , [65, 20]         , ...
                           'PlotBoxAspectRatio', [4.40, 1.86, 1.00], ...
                           'Projection'        , 'Perspective'     , ...
                           'Box'               , 'on');
                       
%% Options
% Mesh
if removeS
    if flatMesh
        opt.mesh = {mesh_size, ['avgThick/G_benchmark_composite_cellSize' ...
                                num2str(mesh_size) '_Sremoved_v2.mat']};
    else
        opt.mesh = {mesh_size, ['G_benchmark_composite_cellSize' ...
                                num2str(mesh_size) '_Sremoved_ThickVar.mat']};
    end
else
    opt.mesh = {mesh_size, ['G_benchmark_composite_cellSize' num2str(mesh_size) '_v2.mat']};
end

%BCs
opt.topboundary = 'curved';  % 'straight' or 'curved'
opt.bctype = 'cp';           % cp: constant pressure at top; pvm: pore volume mult

% Other conditions
opt.temp_C = temp_C;
opt.modelParams = modelParams;

% Injection protocol
if inj_case == 1 || inj_case == 2
    opt.wellno = 4;                                  % n ports total                  
    opt.rate = {37.5, 1;     % {mL/min, times (sec);  inj 1
                0, 1;        %   "   ,   "   ; rampdown 1
                37.5, 1;     %  mL/min, times (sec)}  inj 2
                0, 1;
                37.5, 1;     %  mL/min, times (sec)}  inj 2
                0, 1};       %  seconds, rampdown 2     
    opt.schedule = [0, 30, 60, 90, 120, 150, 180];     % [startI1, endI1, startI2, endI2, end sim]                                      
    opt.hyster = 'off';
    opt.vapoil = false;         % must be set, and must be true or false

elseif inj_case == 3
    opt.injloc = [0.012, 0.93, 1.01;
                  0.013, 1.73, 0.61];
    opt.ploc = [0.012, 1.53, 0.8;
                0.012, 1.73, 0.2];
    opt.wellno = 2;                              % n ports total   
    if strcmp(rampdown, 'fast') 
        % 5.42g total in I1, 2.71g total injected in I2
        %irate1 = 5.42/(5*3600)*1e-3*(1/1.7839549); % m3/s at surface conditions
        %irate2 = 2.71/(2*3600 + 45*60)*1e-3*(1/1.7839549);

        % 5.51g in I1, 3.03 in I2 according to NIST standard
        dCO2 = 1.8683; % 1.013 bar, 15.56C
        if temp_C == 2525 % surface std at 25C
            dCO2 = 1.7839549318418;
        end
        irate1 = 5.5157/(5*3600)*1e-3*(1/dCO2); % m3/s at surface conditions
        irate2 = 3.0337/(2*3600 + 45*60)*1e-3*(1/dCO2);
        mrate = [1e-2 2e-2 5e-2 0.1 0.2 0.5 1];
        opt.rate = {mrate*irate1, 1:7;
                    fliplr([0 mrate(1:6)])*irate1, 1:7;
                    mrate*irate2, 1:7;
                    fliplr([0 mrate(1:6)])*irate2, 1:7};
%         opt.rate = {[0.1 0.2 0.5 1 2 5 10]*minj, 2:2:14;    % {mL/min, s;  inj 1
%                     [5 2 1 0.5 0.2 0.1 0]*minj, 2:2:14;     %   "   ,   "   ; rampdown 1
%                     [0.1 0.2 0.5 1 2 5 10]*minj, 2:2:14;    %  mL/min, times (s)}  inj 2
%                     [5 2 1 0.5 0.2 0.1 0]*minj, 2:2:14};    % rampdown 2
        opt.schedule = [0, 300, 135, 300, 5*24*60];
            
    else
        opt.rate = {[0.1 0.2 0.5 1 2 5 10]*minj, 2:2:14;    % {mL/min, s;  inj 1
                    [6 2 1.5 0.45 0.1 0.01 0]*minj, 2:2:14;     %   "   ,   min   ; rampdown 1
                    [0.1 0.2 0.5 1 2 5 10]*minj, 2:2:14;    %  mL/min, times (s)}  inj 2
                    [6 2 1.5 0.45 0.1 0.01 0]*minj, 2:2:14};    % rampdown 2
        opt.schedule = [0, 298, 135, 298, 5*24*60];
    end
    opt.hyster = 'on';
    opt.vapoil = false;         % must be set, and must be true or false
    opt.D = D;                  % diffusion coeff. [m2/s]
    opt.minSat = sgmin;         % def 0.02
    opt.pcD = pcD;
    opt.dF = dF;
    opt.rampdown = rampdown;
end


if nargin > 12       % for welltest match only
    if strcmp(batch, 'all') || isnumeric(batch)
        nrock = 7;
        nvar = 3;
        nsim = nvar^nrock;
        M = zeros(nsim, nrock);
        M(:,nrock) = repmat([0.8, 1, 1.5]', nsim/nvar, 1);           % G
        M(:,nrock-1) = repmat([repelem(0.7, nvar), repelem(1, nvar), repelem(1.3, nvar)]', ...
            nsim/(nvar^2), 1);                     % F
        M(:,nrock-2) = repmat([repelem(0.7, nvar^2), repelem(1, nvar^2), repelem(1.3, nvar^2)]', ...
            nsim/(nvar^3), 1);                     % E
        M(:,nrock-3) = repmat([repelem(0.7, nvar^3), repelem(1, nvar^3), repelem(1.3, nvar^3)]', ...
            nsim/(nvar^4), 1);                     % D
        M(:,nrock-4) = repmat([repelem(0.6, nvar^4), repelem(1, nvar^4), repelem(1.4, nvar^4)]', ...
            nsim/(nvar^5), 1);                     % C
        M(:,nrock-5) = repmat([repelem(0.5, nvar^5), repelem(1, nvar^5), repelem(1.5, nvar^5)]', ...
            nsim/(nvar^6), 1);                     % ESFsup
        M(:,nrock-6) = repmat([repelem(0.5, nvar^6), repelem(1, nvar^6), repelem(1.5, nvar^6)]', ...
            1, 1);                                 % ESF
        size_batch = nsim/nvar^2;
        id = [0 size_batch*(1:nvar^2)];
        if strcmp(batch, 'all')
            m = M;
        else
            m = M(id(batch)+1:id(batch+1), :);
        end
    
    else    % 2nd matching iteration
        nrock = 4;
        nvar = 4;
        id_var = [1, 3, 4, 7];
        id_set = [2 5 6];
        val_set = [1.5, 0.7 0.7];
        nsim = nvar^nrock;
        M = zeros(nsim, nrock+numel(id_set));
        M(:,id_set) = repmat(val_set, nsim, 1);
        M(:,id_var(end)) = repmat([0.7 0.8, 0.9, 1]', nsim/nvar, 1); 
        M(:,id_var(3)) = repmat([repelem(0.7, nvar), repelem(0.9, nvar), ...
                                 repelem(1.1, nvar),  repelem(1.3, nvar)]', ...
                         nsim/(nvar^2), 1); 
        M(:,id_var(2)) = repmat([repelem(0.6, nvar^2), repelem(0.8, nvar^2), ...
                                 repelem(1.2, nvar^2), repelem(1.4, nvar^2)]', ...
                         nsim/(nvar^3), 1); 
        M(:,id_var(1)) = repmat([repelem(0.5, nvar^3), repelem(0.667, nvar^3), ...
                                 repelem(0.834, nvar^3), repelem(1, nvar^3)]', ...
                         nsim/(nvar^4), 1); 
        size_batch = nsim/2;
        id = [0 size_batch*(1:2)];
        if strcmp(batch, 'all2')
            m = M;
        elseif strcmp(batch, 'first')
            m = M(id(1)+1:id(1+1), :);
        elseif strcmp(batch, 'second')
            m = M(id(2)+1:id(2+1), :);
        end
    end
    
end   

end
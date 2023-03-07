function [mesh, opt, wells, time, resPlots] = simOpts(scenario, studyCase)
%
%
%

% Water column 
opt.watCol   = 50;                                                          % [m] seafloor water column

% Boundary Conditions
opt.bc       = 'pvMult';                                                    % 'noflow' or 'mixed' (p on sides)

% Acceleration
opt.mex      = 1;                                                           % 0 or 1, set up mex-acceleration
opt.compiled = 1;                                                           % 0 or 1, set up compiled ls

% Wells and simulation schedule
wells.num    = 3;                                                           % Number of wells in domain
%time         = [1, 2]*year;                                             

% Results plots default
resPlots    = 'nextToFault';

if scenario == 1    
    if studyCase == 1   
        time         = [20, 100]*year;                                      % [y] injTime, simTime
        % Mesh
        mesh.type    = 'coarse';                                            % 'coarse' or 'fine'
        mesh.reduce  = 0;
        
        % Rock and fluid physics
        opt.petr      = {'isotr-homog','anisotropic', 'sand0_clayLogAvg'};   % {1}: poro, {2}: permd, {3}: Vsh units
        opt.fk_pred   = 'spe02';                                             % fault k predictive model
        opt.fSandPoro = 'constant';                                          % pure sand porosity for mixing model; see endMemberPoro.m
        opt.fluid     = 'krPcBase_pvtBase';                                  % curves and PVT to use
        
        opt.FaultPc  = 'seal';                                              % Fault Pc model. 'seal' is assigned same as for seal units
        
        wells.loc    = [nan, 13200, 2062];                                  % [m] ([x],y,z) location of the wells
        %wells.loc    = [nan, 19688, 2546.8];
        wells.mrate  = 1*10^9/(year*wells.num);                             % 1 Mt/y to [kg/s]
        
        % Results Plots
        resPlots = 'nextToFault';                                           % 'nextToFault' or 'center' or []
        
    elseif studyCase == 2 % base case sc1 (next to fault & center), DEFSIM
        wells.num    = 1;
        opt.studyCase = studyCase;
        time         = [50, 500]*year;                                      % [y] injTime, simTime
        mesh.type    = 'sc2';
        mesh.reduce  = 1; 
        mesh.ampLyr  = 21;
        thi          = [linspace(48.5,400,20) 450 600 750 1000 1600:200:2600 2600]; % 62 + 1 layers total (2 sides)
        mesh.thick   = [fliplr(thi) 30 thi];
        opt.sc       = 1;
        opt.petr     = {'isotr-homog','anisotropic', 'sandLogAvg_clayLogAvg'};
        opt.fk_pred  = 'spe02';
        opt.fSandPoro = 'linearInterp';
        opt.fluid    = 'sc1_krHyst_revised';        
        opt.faultPc  = 'scaled';                                            % based on S, poro and perm
        wells.loc    = [nan, 13200, 2062];        
        wells.mrate   = 1*10^9/(year);                                      % 1 Mt/y to [kg/s] (per well)
        resPlots      = 'nextToFault_sc2';
        %resPlots     = 'center_meshSc2';    
    elseif studyCase == 3
        wells.num    = 1;
        opt.studyCase = studyCase;
        time         = [50, 500]*year;                                      % [y] injTime, simTime
        mesh.type    = 'sc2';
        mesh.reduce  = 1; 
        mesh.ampLyr  = 21;
        thi          = [linspace(48.5,400,20) 450 600 750 1000 1600:200:2600 2600]; % 62 + 1 layers total (2 sides)
        mesh.thick   = [fliplr(thi) 30 thi];
        opt.sc       = 1;
        opt.petr     = {'isotr-homog','anisotropic', 'sandLogAvg_clayLogAvg'};
        opt.fk_pred  = 'spe02';
        opt.fSandPoro = 'linearInterp';
        opt.fluid    = 'sc1_krHyst_revised';        
        opt.faultPc  = 'scaled';                                            % based on S, poro and perm
        wells.loc    = [nan, 12810, 2060];        
        wells.mrate   = 1*10^9/(year);                                      % 1 Mt/y to [kg/s] (per well)
        resPlots      = 'nextToFault_sc2';                                % 'nextToFault' or 'center' or []
    elseif studyCase == 4
        wells.num    = 1;
        opt.studyCase = studyCase;
        time         = [50, 500]*year;                                      % [y] injTime, simTime
        mesh.type    = 'sc2';
        mesh.reduce  = 1; 
        mesh.ampLyr  = 21;
        thi          = [linspace(48.5,400,20) 450 600 750 1000 1600:200:2600 2600]; % 62 + 1 layers total (2 sides)
        mesh.thick   = [fliplr(thi) 30 thi];
        opt.sc       = 1;
        opt.petr     = {'isotr-homog','anisotropic', 'sandLogAvg_clayLogAvg'};
        opt.fk_pred  = 'spe02';
        opt.fSandPoro = 'linearInterp';
        opt.fluid    = 'sc1_krHyst_revised';        
        opt.faultPc  = 'scaled';                                            % based on S, poro and perm
        wells.loc    = [nan, 19900, 2400];        
        wells.mrate   = 1*10^9/(year);                                      % 1 Mt/y to [kg/s] (per well)
        resPlots      = 'center_meshSc2'; 
    end
    
elseif scenario == 2
    if strcmp(studyCase, 'test')
        time          = [30, 200]*year;                                      % [y] injTime, simTime
        mesh.type     = 'sc21';
        mesh.reduce   = 1;
        thi           = [linspace(50,50,40) linspace(100,100,25) ...
                       300 500 700 1100 1600:200:2600 2800];
        mesh.thick    = [fliplr(thi) thi];
        opt.petr      = {'isotr-homog','anisotropic', 'sand0_clayLogAvg'};   % {1}: poro, {2}: permd, {3}: Vsh units
        opt.fk_pred   = 'equiv';                                             % fault k predictive model
        opt.fSandPoro = 'constant';                                          % pure sand porosity for mixing model; see endMemberPoro.m
        opt.fluid     = 'krBase_pcFaultScaled_pvtBase';                      % curves and PVT to use
        opt.faultPc   = 'scaled';                                            % Fault Pc model. 'seal' is assigned same as for seal units
        wells.loc     = [nan, 13200, 2062];                                  % [m] ([x],y,z) location of the wells
        wells.mrate   = 1*10^9/(year*wells.num);                             % 1 Mt/y to [kg/s]
        resPlots      = 'nextToFault';                                       % 'nextToFault' or 'center' or []
    
    elseif strcmp(studyCase, 'sc2_baseCase')  % base case sc2 with spe02 and hyst, DEFSIM
        opt.studyCase = studyCase;
        time          = [50, 500]*year; 
        mesh.type     = 'sc2';
        mesh.reduce   = 1;
        mesh.ampLyr   = 21;
        %thi          = [linspace(50,175,40) 300 500 700 1100 1600:200:2600  2800];   % 102 layers total  (2 sides)
        thi           = [linspace(48.5,400,20) 450 600 750 1000 1600:200:2600 2600]; % 62 + 1 layers total (2 sides)
        mesh.thick    = [fliplr(thi) 30 thi];
        opt.sc        = 2;
        opt.petr      = {'isotr-homog','anisotropic', 'sandLogMin_clayLog70'};  % {1}: poro, {2}: permd, {3}: Vsh units
        opt.fk_pred   = 'spe02';   %smearsUpscaled                            % fault k predictive model
        opt.fSandPoro = 'linearInterp';                                       % pure sand porosity for mixing model; see endMemberPoro.m
        %opt.fluid     = 'krHyst_pcFaultScaled_pvtBase';                       % curves and PVT to use for kr hysteresis
        opt.fluid     = 'sc2_krHyst_revised';                                 % revised, consistent with sc1
        opt.faultPc   = 'scaled';                                             % Fault Pc model. 'seal' is assigned same as for seal units
        opt.SGR       = 'auto';                                               % 'auto' uses faultPermPred.m
        wells.loc     = [nan, 13200, 2062];                                   % [m] ([x],y,z) location of the wells
        wells.mrate   = 1*10^9/(year*wells.num);                              % 1 Mt/y to [kg/s]
        resPlots      = 'nextToFault_sc2';
        
    elseif strcmp(studyCase, 'sc2_predict_theta30_coarse_nohyst')
        opt.studyCase = studyCase;
        time          = [30, 200]*year; 
        mesh.type     = 'sc2';
        mesh.reduce   = 1;
        mesh.ampLyr   = 21;
        thi           = [linspace(50,1000,8) 1500:500:2500 ...
                         3300 4200 4785]; % 28 + 1 layers total (2 sides)
        mesh.thick    = [fliplr(thi) 30 thi];
        opt.sc        = 2;
        opt.petr      = {'isotr-homog','anisotropic', 'sandLogMin_clayLog70'};  % {1}: poro, {2}: permd, {3}: Vsh units
        opt.fk_pred   = 'predict';   %smearsUpscaled                            % fault k predictive model
        opt.fSandPoro = 'predict';                                       % pure sand porosity for mixing model; see endMemberPoro.m
        opt.fluid     = 'noHyst_FaultPcKrUps';                       % curves and PVT to use for kr hysteresis
        opt.faultPc   = 'section_predict';                                             % Fault Pc model. 'seal' is assigned same as for seal units
        wells.loc     = [nan, 13200, 2062];                                   % [m] ([x],y,z) location of the wells
        wells.mrate   = 1*10^9/(year*wells.num);                              % 1 Mt/y to [kg/s]
        resPlots      = 'nextToFault_sc2';
    
    elseif strcmp(studyCase, 'sc2_predict_theta30_coarse')
        opt.studyCase = studyCase;
        time          = [30, 200]*year; 
        mesh.type     = 'sc2';
        mesh.reduce   = 1;
        mesh.ampLyr   = 21;
        thi           = [linspace(50,1000,8) 1500:500:2500 ...
                         3300 4200 4785]; % 28 + 1 layers total (2 sides)
        mesh.thick    = [fliplr(thi) 30 thi];
        opt.sc        = 2;
        opt.petr      = {'isotr-homog','anisotropic', 'sandLogMin_clayLog70'};  % {1}: poro, {2}: permd, {3}: Vsh units
        opt.fk_pred   = 'predict';   %smearsUpscaled                            % fault k predictive model
        opt.fSandPoro = 'predict';                                       % pure sand porosity for mixing model; see endMemberPoro.m
        opt.fluid     = 'krHyst_FaultPcKrUps';                       % curves and PVT to use for kr hysteresis
        opt.faultPc   = 'section_predict';                                             % Fault Pc model. 'seal' is assigned same as for seal units
        wells.loc     = [nan, 13200, 2062];                                   % [m] ([x],y,z) location of the wells
        wells.mrate   = 1*10^9/(year*wells.num);                              % 1 Mt/y to [kg/s]
        resPlots      = 'nextToFault_sc2';
        
    elseif strcmp(studyCase, 'sc2_predict_theta30')
        opt.studyCase = studyCase;
        time          = [30, 200]*year; 
        mesh.type     = 'sc2';
        mesh.reduce   = 1;
        mesh.ampLyr   = 21;
        thi           = [linspace(48.5,400,20) 450 600 750 1000 1600:200:2600 2600]; % 62 + 1 layers total (2 sides)
        mesh.thick    = [fliplr(thi) 30 thi];
        opt.sc        = 2;
        opt.petr      = {'isotr-homog','anisotropic', 'sandLogMin_clayLog70'};  % {1}: poro, {2}: permd, {3}: Vsh units
        opt.fk_pred   = 'predict';   %smearsUpscaled                            % fault k predictive model
        opt.fSandPoro = 'predict';                                       % pure sand porosity for mixing model; see endMemberPoro.m
        opt.fluid     = 'krHyst_FaultPcKrUps';                       % curves and PVT to use for kr hysteresis
        opt.faultPc   = 'section_predict';                                    % Fault Pc model. 'seal' is assigned same as for seal units
        wells.loc     = [nan, 13200, 2062];                                   % [m] ([x],y,z) location of the wells
        wells.mrate   = 1*10^9/(year*wells.num);                              % 1 Mt/y to [kg/s]
        resPlots      = 'nextToFault_sc2';
        
    elseif strcmp(studyCase, 'sc2_predict_theta30_nohyst')
        opt.studyCase = studyCase;
        time          = [30, 200]*year; 
        mesh.type     = 'sc2';
        mesh.reduce   = 1;
        mesh.ampLyr   = 21;
        thi           = [linspace(48.5,400,20) 450 600 750 1000 1600:200:2600 2600]; % 62 + 1 layers total (2 sides)
        mesh.thick    = [fliplr(thi) 30 thi];
        opt.sc        = 2;
        opt.petr      = {'isotr-homog','anisotropic', 'sandLogMin_clayLog70'};  % {1}: poro, {2}: permd, {3}: Vsh units
        opt.fk_pred   = 'predict';   %smearsUpscaled                            % fault k predictive model
        opt.fSandPoro = 'predict';                                       % pure sand porosity for mixing model; see endMemberPoro.m
        opt.fluid     = 'noHyst_FaultPcKrUps';                       % curves and PVT to use for kr hysteresis
        opt.faultPc   = 'section_predict';                                             % Fault Pc model. 'seal' is assigned same as for seal units
        wells.loc     = [nan, 13200, 2062];                                   % [m] ([x],y,z) location of the wells
        wells.mrate   = 1*10^9/(year*wells.num);                              % 1 Mt/y to [kg/s]
        resPlots      = 'nextToFault_sc2';
        
    end
    
    
elseif scenario == 3
    
end
                                           
end
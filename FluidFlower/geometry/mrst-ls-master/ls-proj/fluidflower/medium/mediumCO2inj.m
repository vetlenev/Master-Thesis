%% Perform CO2 injection (blackoil module)
%--------------------------------------------------------------------------
% medium FluidFlower analysis
%--------------------------------------------------------------------------

clear, close all
mrstModule add upr ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils 
mrstVerbose on


%% Options
inj_type =   2;  % 1 (AC02), 2 (AC07) or 3 (BC01)
mesh_size =  4;  % (fine) 1, 2, 4, 8 (coarse); refers to h (mm)
model_case = 3;  % 1 (least data) to 4 (most data)
inj_mult   = [0.85, 0.85];
%inj_mult   = [0.8, 0.7];  % modif inj by that factor with respect to what was desired experimentally.  
D = 1e-9;
%Fpc = 0;
%Epc = 0.
%folderName = 'exp1_mesh4_modelcase1_D3e-9_I1m0.8_I2m0.7_FPc0_EPc0.1_ESF0.1mm_ESFPce1e-5_F2.8mm_C0.9mm_E1.7mm'; 
if model_case == 1
    m = [];
    folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
                 '_modelcase' num2str(model_case) '_D' num2str(D)...
                 '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
                 '_ESF0.2mm_ESFPceSg1e-5_C_0.66mm_E_1.45mm_F_1.77mm']; % 1st run
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                  '_modelcase' num2str(model_case) '_D' num2str(D)...
%                  '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                  '_ESF0.2mm_ESFPceSg1e-4_C_0.66mm_E_1.45mm_F_1.77mm' ...
%                  '_FPc0_EPc0.125_CPc0.5'];
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                  '_modelcase' num2str(model_case) '_D' num2str(D)...
%                  '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                  '_ESF0.1mm_ESFPceSg1e-4_C_0.66mm_Cf_0.2mm_E_1.45mm_Fsupinf_3mm' ...
%                  '_Fmed_2mm_FPc0_EPc0.125_Cpc0.5'];
    topDir = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\medium\';
elseif model_case == 2
    m = [1 1 1 1 1 1 1];                     % k mult initial
    %m = [1 1 1/3 1.5 1.5 1.5 1];            % k mult [ESF, C, Cf, E, Fsup, Finf, Fmid]
    folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
                  '_modelcase' num2str(model_case) '_D' num2str(D) ...
                  '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
                  '_repPhi_repk_repPc_model1krgkrw']; % initial run
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_model1krgkrw_FPc0'];
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_model1krgkrw_FPc0_Epc0.25_Cpc0.5'];
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_model1krgkrw_FPc0_Epc0.125_Cpc0.33'];
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_model1krgkrw_FPc0_Epc0.18_Cpc0.33_Fk2'];
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc0.33_Fksupinf1.5'...
%                   '_Fkmid1_Ek1.5_Cfk0.33'];
    topDir = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\medium\model2\';
elseif model_case == 3
    %m = [1 1 1 1 1 1 1];                % initial
    m = [1/3 1 1/4 1.2 1.7 1.7 1.1];    % match for D=1e-9
    %m = [1/3 1 1/4 1.6 2.3 2.3 1];       % match for D=3e-9
    %m = [1/3 1 1/4 1.8 2.5 2.5 1];       % match for D=5e-9
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_repPc_repSwcSgt_model1krgkrw']; % initial run
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_repSwcSgt_model1krgkrw_Cpc1.5'];
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33_Fk2'];
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%                   '_Fk1.5_Cfk0.5'];
    folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
                  '_modelcase' num2str(model_case) '_D' num2str(D) ...
                  '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
                  '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
                  '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2']; % match for D=1e-9
%     folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%                   '_Fk2.3_Fkmid1_Cfk0.25_Ek1.6']; % match for D=3e-9
%         folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%                   '_Fk2.5_Fkmid1_Cfk0.25_Ek1.8']; % match for D=5e-9
    
    topDir = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\medium\model3\';
    %topDir = '/home/lsalo/matlab/sim_data/mrst/fluidflower/medium/AC02/';
end
if inj_type == 2
    topDir = '/home/lsalo/matlab/sim_data/mrst/fluidflower/medium/AC07/';
    if model_case == 1
%         folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                  '_modelcase' num2str(model_case) '_D' num2str(D)...
%                  '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                  '_ESF0.1mm_ESFPceSg1e-4_C_0.66mm_Cf_0.2mm_E_1.45mm_Fsupinf_3mm' ...
%                  '_Fmed_2mm_FPc0_EPc0.125_Cpc0.5'];
%         folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%                  '_modelcase' num2str(model_case) '_D' num2str(D)...
%                  '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                  '_ESF0.1mm_ESFPceSg1e-4_C_0.66mm_Cf_0.2mm_E_1.45mm_Fsupinf_3mm' ...
%                  '_Fmed_2mm_FPc0_EPc0.125_Cpc0.5_CfinfPce7.3mb_CfsupPce4.5mb'];
        folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
                 '_modelcase' num2str(model_case) '_D' num2str(D)...
                 '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
                 '_ESF0.1mm_ESFPceSg1e-4_C_0.66mm_Cf_0.2mm_E_1.45mm_Fsupinf_3mm' ...
                 '_Fmed_2mm_FPc0_EPc0.125_Cpc0.5_CfinfPce5mb_CfsupPce3.5mb'];
    elseif model_case == 2
%         folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc0.33_Fksupinf1.5'...
%                   '_Fkmid1_Ek1.5_Cfk0.33'];
%           folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%                   '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                   '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                   '_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc0.33_Fksupinf1.5'...
%                   '_Fkmid1_Ek1.5_Cfk0.33_CfinfPce7.3mb_CfsupPce4.5mb'];
          folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
                  '_modelcase' num2str(model_case) '_D' num2str(D) ...
                  '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
                  '_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc0.33_Fksupinf1.5'...
                  '_Fkmid1_Ek1.5_Cfk0.33_CfinfPce5mb_CfsupPce3.5mb'];
    elseif model_case == 3
%         folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%             '_modelcase' num2str(model_case) '_D' num2str(D) ...
%             '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%             '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%             '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2'];
%         folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%             '_modelcase' num2str(model_case) '_D' num2str(D) ...
%             '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%             '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%             '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce7.3mb_CfsupPce4.5mb'];
%         folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%             '_modelcase' num2str(model_case) '_D' num2str(D) ...
%             '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%             '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%             '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce4.5mb_CfsupPce4.5mb'];
%         folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%             '_modelcase' num2str(model_case) '_D' num2str(D) ...
%             '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%             '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%             '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce4.5mb_CfsupPce4.5mbSge1e-4'];
%         folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%             '_modelcase' num2str(model_case) '_D' num2str(D) ...
%             '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%             '_repPhi_repSwcSgt_model1krgkrw_ESFpc2_ESFk0.33' ...
%             '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce5mb_CrestPce3.5mbSge1e-4'];
%         folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%             '_modelcase' num2str(model_case) '_D' num2str(D) ...
%             '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%             '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%             '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce5mb_CfsupPce3.5mb'];
        folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
            '_modelcase' num2str(model_case) '_D' num2str(D) ...
            '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
            '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5Sge1e-4_ESFpc2_ESFk0.33' ...
            '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce5mbSge1e-3_CfsupPce3.5mbSge1e-3'];
%         folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%             '_modelcase' num2str(model_case) '_D' num2str(D) ...
%             '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%             '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5Sge1e-5_ESFpc2_ESFk0.33' ...
%             '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce5mbSge1e-3_CfsupPce3.5mbSge1e-4'];
%         folderName = ['exp' num2str(inj_type) '_AC07mesh' num2str(mesh_size) ...
%             '_modelcase' num2str(model_case) '_D' num2str(D) ...
%             '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%             '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5Sge1e-5_ESFpc2_ESFk0.33' ...
%             '_Fk1.7_Fkmid1.1_Cfk0.25_Ek1.2_CfinfPce5mbSge1e-4_CfsupPce3.5mbSge1e-4'];
    end
elseif inj_type == 3
    topDir = 'C:\Users\lsalo\matlab\sim_data\mrst\fluidflower\medium\BC01\';
    %topDir = 'C:\Users\Lluis\matlab\sim_data\mrst\fluidflower\medium\BC01\';
    if model_case == 1
%         folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                     '_modelcase' num2str(model_case) '_D' num2str(D)...
%                     '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                     '_ESF0.1mm_ESFPceSg1e-4_C0.66mm_E1.45mm_F3mm' ...
%                     '_FPc0_EPc0.125_Cpc0.5'];
%          folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                     '_modelcase' num2str(model_case) '_D' num2str(D)...
%                     '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                     '_ESF0.1mm_ESFPceSg1e-4_C0.66mm_E1.45mm_F3mm' ...
%                     '_FPc0_EPc0.125'];
            folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
                '_modelcase' num2str(model_case) '_D' num2str(D)...
                '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
                '_ESF0.1mm_ESFPceSg1e-4_C0.66mm_E1.45mm_F3mm' ...
                '_FPc0_EPc0.125_CPceSg1e-4'];
%             folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                 '_modelcase' num2str(model_case) '_D' num2str(D)...
%                 '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                 '_ESF0.1mm_ESFPceSg1e-4_C0.66mm_E1.45mm_F3mm' ...
%                 '_FPc0_EPc0.125_CPc2_CPceSg1e-5'];
%             folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                 '_modelcase' num2str(model_case) '_D' num2str(D)...
%                 '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                 '_ESF0.1mm_ESFPceSg1e-4_C0.66mm_E1.45mm_F3mm' ...
%                 '_FPc0_EPc0.125_CPc2_CPceSg1e-5_nohyst'];
%             folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                 '_modelcase' num2str(model_case) '_D' num2str(D)...
%                 '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                 '_ESF0.1mm_ESFPceSg1e-4_C0.66mm_E1.45mm_F3mm' ...
%                 '_FPc0_EPc0.125_CPc2_CPceSg1e-4'];

    elseif model_case == 2
        m = [1 1 1 1 1.5 1.5 1.6 1.6];      % match for D=1e-9 
%         folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                       '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                       '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                       '_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc0.33' ...
%                       '_ESFk1_Ck1_Ek1.5_Fk1.6'];
        folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
                      '_modelcase' num2str(model_case) '_D' num2str(D) ...
                      '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
                      '_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc1.3' ...
                      '_ESFk1_Ck1_Ek1.5_Fk1.6'];
% folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                       '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                       '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                       '_repPhi_repk_model1krgkrw_FPc0_Epc0.15_Cpc2' ...
%                       '_ESFk1_Ck1_Ek1.5_Fk1.6'];
    elseif model_case == 3
        m = [1/3 1 1 1 1.2 1.2 1.7 1.7];    % match for D=1e-9 
%         folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
%                       '_modelcase' num2str(model_case) '_D' num2str(D) ...
%                       '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
%                       '_repPhi_repSwcSgt_model1krgkrw_Cpc1.5_ESFpc2_ESFk0.33' ...
%                       '_Ck1_Ek1.2_Fk1.7']; % match for D=1e-9
        folderName = ['exp' num2str(inj_type) '_mesh' num2str(mesh_size) ...
                      '_modelcase' num2str(model_case) '_D' num2str(D) ...
                      '_I1m' num2str(inj_mult(1)) '_I2m' num2str(inj_mult(2))...
                      '_repPhi_repSwcSgt_model1krgkrw_Cpc3.2_ESFpc2_ESFk0.33' ...
                      '_Ck1_Ek1.2_Fk1.7'];
    end
end
%topDir = 'C:\Users\Lluis\matlab\sim_data\mrst\fluidflower\medium\model2\';
%topDir = '/Users/lluis/Documents/MATLAB/sim_data/mrst/fluidflower/medium/';
[plts, opt] = getSimOpts(inj_type, mesh_size, model_case, inj_mult, D);


%% Mesh
[G, G_dat, unit, wellId, meshpt] = getMesh(inj_type, opt, plts, 0);


%% Rock
rock = getRock(inj_type, model_case, G, G_dat, unit, m, plts, 0, 0);


%% Fluids
[fluid, rock] = getFluid(inj_type, model_case, G, G_dat, rock, unit, opt, plts, 0);


%% Initialize
% theta = 0;      
% R = makehgtform('xrotate',pi*theta/180);    % theta deg. about the x axis
gravity reset on
% gravity( R(1:3,1:3)*gravity().' );
[state0, p_r] = initialize(inj_type, G, fluid, plts, 0);


%% Wells
[W, timesteps, wellInx] = getWells(G, G_dat, rock, opt, model_case, inj_type);


%% Model and acceleration
[model, nls] = getModel(G, rock, fluid, state0, opt); 


%% BCs
[bc] = getBC(G, fluid, p_r, opt);


%% Schedule
schedule = getSchedule(inj_type, timesteps, W, bc, opt);


%% Simulation
N = 2;
maxNumCompThreads(N);
nls.LinearSolver.amgcl_setup.nthreads = N;                                  % Specify threads manually
%[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
%                                           'NonLinearSolver', nls, ...
%                                           'Verbose', true);
 
problem = packSimulationProblem(state0, model, schedule, folderName, 'Name', folderName, ...
                                'Directory', topDir, 'NonLinearSolver', nls);
[ok, status] = simulatePackedProblem(problem);
[wellSols, states, report] = getPackedSimulatorOutput(problem);  


%% Results
% CO2 mass and V end of I1 (not accounting for rock porevolume changes)

% Get areas
c_bound = 0.15*1.4; % based on what we can see in plot & pictures
sg_bound = 1e-3;
A_simulation_cm2 = getMediumAreas(G, G_dat, unit, inj_type, states, model, ...
                                  timesteps, sg_bound, c_bound, problem.Name, 1);

% if numel(states) == numel(timesteps)
%     c_bound = [0.05 0.05];
%     t_fingers_s = getFingerTimes(G, G_dat, model, states, unit, ...
%                                  timesteps, c_bound, problem.Name, 1);
% end

% Mass
%[total_mass] = getMass(model, states, timesteps, problem.Name, 1);

% Results plots
sg_bound = 1e-3;
c_bound = 0.15*1.4;
if inj_type == 1
    %       54', 154', 5h15', 7h
    tid = [54, 154, 186, 207];
    %      3h   6h   9h   12h  18h  24h
    %tid = [159, 195, 231, 250, 286, 322];
    plotsMedium(G, rock, model, meshpt, plts, states, sg_bound, ...
                c_bound, tid, timesteps)
elseif inj_type == 2
    tid = 290;
    plotsMedium(G, rock, model, meshpt, plts, states, sg_bound, ...
                c_bound, tid, timesteps)
elseif inj_type == 3
    %plotsBilbo(G, rock, meshpt, plts, s, componentPhaseMass, c_bound, 200, 205)
    %      3h   6h   9h   12h  18h  24h
    tid = [159, 195, 231, 250, 286, 322];
    plotsMedium(G, rock, model, meshpt, plts, states, sg_bound, ...
                c_bound, tid, timesteps)
end

% compute quantitites
model = model.validateModel();
bhp = zeros(1, numel(states));
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - state0.pressure;
    rho = model.getProp(states{n}, 'ComponentPhaseDensity');
    states{n}.FlowProps.ComponentPhaseDensity = rho{:,1}; % 2 components in brine phase.
    %pc = model.getProp(states{n}, 'CapillaryPressure');
    %states{n}.FlowProps.CapillaryPressure = pc{1,2};
    %states{n}.FlowProps.RelativePermeability = model.getProp(states{n}, 'RelativePermeability');
    bhp(n) = wellSols{n}.bhp;
end
states{1}.reg = model.rock.regions.saturation;

% basic overview
%rescells = find(ismember(G_dat.p, unit.F));
%ESFcells = find(ismember(G_dat.p, unit.ESF));
plts.fig3D(); plotToolbar(G, states, 'edgealpha', 0.2); hold on
%plotGrid(G, rescells, 'edgecolor', [0.7 0.7 0.7]);
%plotGrid(G, ESFcells, 'edgecolor', [0.4 0.4 0.4]);
%plot(meshpt.boundary(:,1), meshpt.boundary(:,2), 'k')
xmx = max(G.faces.centroids(:,1));
for n=1:numel(meshpt.lines)
    if n < 5
        clr = [0.3 0.3 0.3];
    else
        clr = [0.7 0.7 0.7];
    end
    xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
    plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
end
%plotLinePath(meshpt.lines(5:end), 'b');
%plotLinePath(meshpt.lines(1:4), 'r');
%cmap = getMyCmap('seashore');
cmap = flipud(cmocean('deep'));
plts.setAxProps(gca), colormap(cmap), c = colorbar; caxis([1e-6 1])
axis equal off
view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
set(gca, 'ColorScale', 'log')
%c.Limits = [1e-6 1];

% video
dir = 'C:\Users\lsalo\matlab\';
type = 'cgw';
idt = [5:5:150 154 155:300];
%idt = [5:5:150 154];
getStatesVideo(dir, type, folderName, idt, plts, ...
               G, rock, states, model, timesteps, meshpt);

% % well cell pressure
% dpw = zeros(numel(states), numel(wellId));
% for n=1:numel(states)
%     dpw(n, :) = states{n}.pressure(wellId) - state0.pressure(wellId);
% end
% figure(2)
% plot(cumsum(timesteps(1:numel(states)))/minute, dpw)
% legend('P1', 'P2', 'I1', 'P3', 'I2', 'P4', 'P5', 'P6')
% grid on

% Plot cell data at a given time
% idt = 88;
% plts.fig3D();
% plotCellData(G, states{idt}.s(:,2), 'edgecolor', 'none')
% %plotCellData(G, states{idt}.PVTProps.Density{1}, 'edgecolor', 'none')
% hold on
% xmx = max(G.faces.centroids(:,1));
% for n=1:numel(meshpt.lines)
%     if n < 5
%         clr = [0.3 0.3 0.3];
%     else
%         clr = [0.7 0.7 0.7];
%     end
%     xcrd = repelem(xmx, size(meshpt.lines{n}, 1));
%     plot3(xcrd', meshpt.lines{n}(:,1), meshpt.lines{n}(:,2), '-', 'color', clr)
% end
% plts.setAxProps(gca)
% colormap(hot)
% %colormap(cmocean('balance'))
% c = colorbar; 
% ylabel(c, '$S_{g}$ [-]', 'fontsize', 14, 'interpreter', 'latex')
% %ylabel(c, '$\rho_\mathrm{aq}$ [kg/m$^3$]', 'fontSize', 14, 'interpreter', 'latex')
% %caxis([995.1 995.4])
% axis equal off
% view([90 0]), hold off %ylim([0.42 0.48]), zlim([0.40 0.47])
% set(gca, 'ColorScale', 'log')
% caxis([1e-4 1])

% Brine Density at a given time
% id = states{end}.rs > 0.001;
% plotGrid(G,'edgealpha', 0.2)
% plotCellData(G, states{end}.PVTProps.Density{1}(id), find(id), 'edgealpha', 0.2)
% plts.setAxProps(gca), colormap(turbo), c = colorbar;
% axis equal off, view([90 0])

cmap = cmap/255;
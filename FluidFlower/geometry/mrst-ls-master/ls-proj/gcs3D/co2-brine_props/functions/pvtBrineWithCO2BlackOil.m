function [t, rho_co2_s, rho_brine_s] = pvtBrineWithCO2BlackOil(T, P, S, ...
                                        saltVar, iterate, figs, directory)
%
% Lluís Saló-Salgado (lsalo@mit.edu), September 2019
% Updates:
% LS, July 2022:  adjust max p vals and update plots
%
% DESCRIPTION:
% Compute and save a table with R_s, P_aq, B_aq, mu_aq for a saline aqueous
% solution with dissolved CO2, and another table with (R_s), P_gas, B_gas 
% and mu_gas for a CO2-rich gaseous phase (with or without dissolved H2O). 
% Intended for usage with a black-oil formulation. 
%
% MAIN REFERENCE:
% Hassanzadeh et al., IJGGC (2008). There are a few typos in the paper,
% look for "typo" in the comments below.
%
% INPUT:
% T:    cell indicating reference temperature (typically temperature at
%       reservoir depth, since it is where it matters most. Note that the 
%       black-oil formulation is isothermal) and units. Accepted C, K, F.
%       Structured as a {'key',value} pair. 
%       Example: T = {'C', 100};
%
% P:    cell indicating pressure values at which Rs and Bw should be
%       computed, units (Accepted Pa, MPa, bar, psi), and type.
%       type: 'mMn'  --> [min p, max p, num vals]
%             'vals' --> array of size (1, n) containing all pressure
%                        values of interest.
%       Structured as {'units','type', array}. 
%       Example: pval = {'Pa','mMn',[10^5, 10^8, 12]};
%
% S:    cell indicating reference salinity and units (accepted 'ppm' 
%       (mg salt/kg solution), 'g/kg' (g salt /kg solution), 'sal' (ppt in 
%       g salt / kg solution), 'ppml' (mg salt / l solution)).
%       Structured as {'units','type', value}. 'type' indicates the type of
%       fluid for ion and salt breakdown assignment, based on the input
%       total salinity value. Accepted options for 'type' are 'NaCl',
%       'seawater' and 'brine' for nonmolal 'unit' and 'NaCl' or 'CaCl2'
%       for 'm' unit.
%       Example: S = {'ppm', 'brine', 10^5};
%
% saltVar: double indicating whether the change in salt mole fraction due
%          to CO2 dissolution into the brine should be accounted for.
%          Example: saltVar = 0; (0 for "no" or 1 for "yes")
%
% iterate: double indicating whether iteration should be performed to
%          account for water vaporization in the gasous (CO2 rich) phase,
%          when determining the interaction parameters of H2O and CO2 for
%          computing points 1 to 6.
%          Since Y_h2o is usually very small, this shouldn't lead to major
%          changes in usual CO2 storage settings.
%          Example: iterate = 0; (0 for "no" or 1 for "yes")
%
% figs:    double indicating whether figures illustrating the results
%          should or should not be plotted.
%          Example: figs = 0; (0 for "no" or 1 for "yes")
%
% Optional argument:
% directory:  string indicating the full path to the directory where output
%             table file should be saved.
%             Example: 'C:\Users\lsalo\matlab\mrst-2019a\myutils\'
%
% OUTPUT:
% t: Structure specifying temperature, pressures, brine salinity and
%    aqueous and gase phases values of interest. These are (see below for
%    naming conventions):
%    t.aq:  m_co2, x_co2, X_co2, X_salt, rho_b, rho_aq, Rs, B and mu
%    t.gas: Z_co2, phi_co2, gamma_co2, y_h2o, Y_h2o, rho_g, Rv, B and mu
%
% If directory is passed:
% Saves a table in .txt format in the specified directory. Units are 
% ECLIPSE 'METRIC' units. The table is saved as:
% - Aqueous phase: a table called fileName of size (n, 4), where n is the 
%   number of pressure values of interest. Rs is the 1st column, p is the 
%   2nd column and the corresponding B_aq and viscosity are the 3rd and 4th
%   columns. fileName = PVTO_Tval_PminPMax_StypeValUnit.txt
% - Gaseous phase: a table called fileName of size (n, 3/4) with the same
%   format as above for the aqueous phase. In case Rv = 0 (dry gas),
%   fileName = PVDG_... and size is (n, 3); in case Rv \neq 0 (wet gas) 
%   fileName = PVTG_... and size is (n, 4).
%
% NOTES:
% - Two additional values of pressure above maximum input pressure are
%   added in order to provide two values for undersaturated aqueous
%   (gaseous) phase(s) for the highest Rs (Rv). This is done because it is
%   the standard ECLIPSE format for .DATA files. For these two values, the
%   gas compressibility factor and fugacity coefficients are not computed
%   since they are not needed, and that is why they will appear as NaN in
%   the output table file.
%
%
% VARIABLE NAME CONVENTIONS:
% T            temperature
% _K _C        Kelvin, Celsius
% P, p         pressure
% S            salinity
% rho          mass density
% rhox         molar density
% mu           viscosity
% Rs, Rv       solution CO2-brine ratio and vaporized H2O-CO2 ratio
% B            formation volume factor
% Z            gas compressibility factor
% phi          fugacity coefficient
% gamma        activity coefficient
% Mw_k         molar mass of k
% m_k          molality of species k in solution (mol k / kg solvent)
% w_k          mass fraction of species k in aqueous phase
% v_k          mass fraction of species k in gaseous phase
% x_k          mole fraction of species k in water phase
% X_k          mole fraction of species k in aqueous phase
% y_k, Y_k     mole fraction of species k in gaseous (CO2 rich) phase
% _b, _brine   brine (H2O + salt(s))
% _gas         gaseous phase (CO2-rich)
% _aq          aqueous saline solution with CO2 (H2O + salt(s) + CO2)
% _s           ECLIPSE standard conditions of T, P (15.56C, 1atm)
% _1, _2       solute and solvent, respectively
% _pt          partial (for partial molar volume)
% nu           salt stoichiometric number
% 
%                       ----------------------

%% EoS
% Pick co2, brine, aqueous (brine + CO2) and gas (CO2 + H2O) EoS for 
% density and viscosity. 
s_unit  = S{1};  s_type = S{2};   s_val  = S{3};

% Density
switch s_type
    case {'NaCl'}
        rho_b_fcn = @(T, P, S) pvtBrineRoweChou1970(T, P, S);             % [kg/m^3]
    otherwise
        rho_b_fcn = @(T, P, S) pvtBrineBW1992(T, P, S);                   % [ " ]
end
rho_co2_fcn = @(T, P, varargin) pvtCO2RK1949(T, P, varargin);             % [ " ]
rho_aq_fcn  = @(Mw_1, Mw_2, x_1, x_2, V_pt, rho_1) ...                    % Garcia (2001)
              (1 + (Mw_1/Mw_2)*(x_1/x_2)) / ((V_pt/Mw_2)*(x_1/x_2) + (1/rho_1));
rho_gas_fcn = @(Mw_gas, V) Mw_gas/V;
  
    
% Viscosity
mu_co2_fcn  = @(T, rho) viscCO2F1998(T, rho);                             % [Pa*s]
mu_aq_fcn   = @(T, P, m_S, w_co2) viscBrineCO2IC2012(T, P, m_S, w_co2);   % [ " ]
mu_gas_fcn  = @(x, Mw, mu) viscGasMixtD1993(x, Mw, mu);                    


%% Check and organize inputs
% ------------------------- Temperature -----------------------------------
T_unit = T{1}; T_val  = T{2};
if strcmp(T_unit, 'K'),     T_K = T_val;
elseif strcmp(T_unit, 'C'), T_K = T_val + 273.15;
elseif strcmp(T_unit, 'F'), T_K = (T_val + 459.67)*(5/9);
else, error("Temperature units must be 'K', 'C' or 'F'.")
end
T_C   = T_K - 273.15;                                                     
T_K_s = 273.15 + 15.56;                                                   

% ---------------------------- Pressure -----------------------------------
p_unit = P{1}; p_val  = P{3}; p_typ  = P{2};
assert(max(strcmp(p_unit, {'Pa','MPa','bar','psi'})))
switch p_typ
    case 'mMn', assert(numel(p_val)==3)
        p = logspace(log10(p_val(1)), log10(p_val(2)), p_val(3))';
    case 'vals'
        p = p_val';
    otherwise, error("p{2} ('type') must be 'mMn' or 'vals'.")
end
if strcmp(p_unit, 'Pa'),        p = p*10^(-5);
elseif strcmp(p_unit, 'MPa'),   p = p*10;
elseif strcmp(p_unit, 'psi'),   p = p*0.06894757;
end
P_bar_s = 1.013529;                                                       
pmax    = max(p);
p       = [p; 1.1*pmax; 1.25*pmax];                                         % for undersaturated oil, highest Rs

% ---------------------------- Salinity -----------------------------------
% Check input 
if numel(S)~=3
        error("number of elements in S cell must be 3")
end

% Initial variables
maxits      = 8;                                                          % max N iterations
saltSpecies = {'NaCl', 'KCl', 'CaSO4', 'CaCl2', 'MgSO4', 'MgCl2'};
%         [ NaCl,     KCl,    CaSO4,  CaCl2,  MgSO4,  MgCl2]
Mw_salt = [58.44277 74.5513 136.1406 110.984 120.3676 95.211]/10^3;       % [kg/mol]
%         [ Na(+),   K(+),  Ca(2+), Mg(2+), Cl(-), SO4(2-)]
Mw_io   = [22.98977 39.0983 40.078 24.305 35.453 96.0626]/10^3;           % ["]

% Compute salt molalities
if ~strcmp(S{1}, 'm')
    switch s_unit
        case {'ppm','ppml'}
            kg_salt = s_val/10^6;                                         % [kg salt / kg solution]
        case {'g/kg', 'sal'}
            kg_salt   = s_val/10^3;                                       % [kg salt / kg solution]                 
        otherwise 
            error("Salinity units must be 'ppm(l)', 'g/l', 'sal', 'm'")
    end
    kg_wat  = 1-kg_salt;                                                  % [kg] water / kg solution
    if strcmp(s_unit,'ppml')                                              % correct for volume
        dif = 1;    tols = 1e-4;    it = 1;
        while dif > tols && it < maxits                                   % find kg water
            rho_bs  = rho_b_fcn(T_K_s, P_bar_s, kg_salt/(kg_wat+kg_salt));
            kg1_w = rho_bs/10^3 - kg_salt;                                % kg h2o / l sol
            dif = abs(kg1_w - kg_wat);
            kg_wat = kg1_w;
            if dif < tols
                break
            elseif it < maxits
                it = it + 1;
                %disp(dif)
            else, error('Iteration to correct for volume did not converge')
            end
        end
    end
    switch s_type
        case 'NaCl'
            m_NaCl  = kg_salt/(kg_wat*Mw_salt(1));                        % [m] 
            m_salt  = [m_NaCl 0 0 0 0 0];                                 % [m]
        case 'seawater'
            % https://www.britannica.com/science/seawater for data
            w_io = [10.68 0.40 0.41 1.28 19.16 2.68]/34.61;
            kg_io = w_io*kg_salt;
            m_io = (kg_io./(Mw_io*kg_wat));                               % [m]
            m_salt = [m_io(1)  m_io(2)  m_io(3)/2  m_io(3)/2 ...
                      m_io(6)-(m_io(3)/2)  m_io(4)-m_io(6)-(m_io(3)/2)]; 
        case 'brine'
            % Data from Hyeong & Capuano (2001);
            w_io = [23.5 0.181 0.676 0.184 36.1 0]/10^3;                  % [kg io / l sol]
            w_io = w_io/sum(w_io);                                        % [-]                     
            kg_io = w_io*kg_salt;   
            m_io = (kg_io./(Mw_io*kg_wat));                               % [m]
            m_salt = [m_io(5)-m_io(2)-2*(sum(m_io(3:4)))  m_io(2)  0  m_io(3) 0  m_io(4)];
        otherwise, error('This salt input is not supported.')
    end
else
    switch s_type
        case 'NaCl'
        m_salt  = [S{3}, 0, 0, 0, 0, 0];                                  % [m]
        kg_salt = S{3}*Mw_salt(1);
        case 'CaCl2'
        m_salt  = [0, 0, 0, S{3}, 0, 0];
        kg_salt = S{3}*Mw_salt(4);
        otherwise, error('This salt is not supported for molal inputs')
    end
    kg_wat = 1;
end

% ------------------- Ion molalities and salt mass fraction -------------------
nu      = [  2,     2,     2,       3,     2,      3];                    
w_salt = kg_salt/(kg_wat+kg_salt);                                        % [-] 
if ~strcmp(s_type,'seawater') && ~strcmp(s_type,'brine')
    m_io   = [m_salt(1) m_salt(2) m_salt(3)+m_salt(4) m_salt(5)+m_salt(6) ...
              sum(m_salt(1:2))+2*(m_salt(4)+m_salt(6)) m_salt(3)+m_salt(5)];
end
   

%% Preliminary variables 
R       = 83.1447;                                          % bar*cm^3/mol*K
Mw_h2o  = 18.01528/10^3;                                    % [kg/mol]
Mw_co2  = 44.0095/10^3;                                     %   "

% Compute mole fractions
m_h2o   = 1/Mw_h2o;                                         % mol H20/ kg H2O
x_salt  = m_salt./(sum(m_salt) + m_h2o);                    % [-] mole fractions (total)
x_h2o   = 1 - sum(x_salt);                                  % [-]
x_salts = m_salt./sum(m_salt);                              % [-] (of salts)
Mw_salts = sum(Mw_salt.*x_salts);                           % [kg/mol] avg

% Compute mixture Mw (Mw_brine) in kg/mol
Mw_b = sum(Mw_salt.*x_salt) + Mw_h2o*x_h2o;                 % [kg/mol]

% Water and CO2 interaction parameters after Spycher et al. (2003)
a_co2       = 7.54*10^7 - 4.13*10^4*T_K;                    % bar*cm^6*K^0.5/mol^2
a_h2o_co2   = 7.89*10^7;                                    % "
b_co2       = 27.8;                                         % cm^3/mol
b_h2o       = 18.18;                                        % "

% Average partial molar volumes of each pure condensed component, they are
% assumed constant in the p interval of interest (Spycher et al., 2003)
Vp_co2      = 32.6;                                         % cm^3/mol (both (g) and (l))
Vp_h2o      = 18.1;                                         % "

% True equilibrium ctnts. at p0=1bar (/kappa params; Spycher et al., 2003)
p0          = 1;                                            % bar
kap0_co2_g  = 10^(1.189 + 1.304*10^(-2)*T_C - 5.446*10^(-5)*T_C^2);
%kap0_co2_l  = 10^(1.169 + 1.368*10^(-2)*T_C - 5.380*10^(-5)*T_C^2);
kap0_h2o    = 10^(-2.209 + 3.097*10^(-2)*T_C - 1.098*10^(-4)*T_C^2 ...
                  + 2.048*10^(-7)*T_C^3);

% Partial molar volume of CO2 (Garcia, 2001)
V_phi     = 37.51 - 9.585*10^(-2)*T_C + 8.740*10^(-4)*T_C^2 ...
            - 5.044*10^(-7)*T_C^3;                          % [cm^3/mol]
V_phi     = V_phi/10^6;                                     % [m^3/mol]

% CO2 and brine molar densities at standard conditions
[~, rhox_co2_s, rho_co2_s] = rho_co2_fcn(T_K_s, P_bar_s);   % [mol/m^3]
% With CoolProp library (EoS from Span & Wagner, 1996)
% rhox_co2_s   = py.CoolProp.CoolProp.PropsSI('DMOLAR',...
%                'T', T_K_s, 'P', P_Pa_s, 'CO2');           % [mol/m^3]
rho_brine_s   = rho_b_fcn(T_K_s, P_bar_s, w_salt);
rhox_brine_s  = rho_brine_s/Mw_b;                           % [mol/m^3]
rho_h2o_s     = pvtBrineRoweChou1970(T_K_s, P_bar_s, 0);    
rhox_h2o_s    = rho_h2o_s/Mw_h2o;                           % ["]


%% Calculate
% Initialize mole fractions in gaseous phase
Y0_h2o      = 0;
Y0_co2      = 1;

% Initialize result table
hdrs.gas    = {'Z_co2', 'phi_co2', 'gamma_co2', 'y_h2o', 'Y_h2o', ...
               'rho_g_kgm3', 'Rv_Sm3Sm3', 'B_g_m3Sm3', 'mu_g_cP'};
hdrs.aq     = {'m_co2', 'x_co2', 'X_co2', 'X_salt', 'rho_b_kgm3', ...
               'rho_aq_kgm3', 'Rs_Sm3Sm3', 'B_b_m3Sm3', 'B_aq_m3Sm3', ...
               'mu_aq_cP'};
vartyp.gas  = cellstr(repmat('double',numel(hdrs.gas),1));
vartyp.aq   = cellstr(repmat('double',numel(hdrs.aq),1));
t.gas       = table('Size',[numel(p), numel(hdrs.gas)], ...
                    'VariableTypes', vartyp.gas, 'VariableNames', hdrs.gas);
t.aq        = table('Size',[numel(p), numel(hdrs.aq)], ...
                    'VariableTypes', vartyp.aq, 'VariableNames', hdrs.aq);
t.P_bar     = p; 
t.T_Kelvin  = T_K; 
t.S_molal.salts = saltSpecies; t.S_molal.vals = m_salt;

% Compute
for n=1:numel(p)
    it = 1;
    P  = p(n);
    if P <= pmax
        while it <= maxits
            % Mixing rules (Prausnitz et al., 1986)
            a_m = Y0_co2^2*a_co2 + Y0_co2*Y0_h2o*a_h2o_co2 + ...
                  Y0_h2o^2*a_h2o_co2;
            b_m = Y0_co2*b_co2 + Y0_h2o*b_h2o;
            
            % 1. Gas molar volume (Redlich and Kwong (1949) EoS, = RK EoS)
            [V, ~, rho_co2] = rho_co2_fcn(T_K, P, a_m, b_m);              % [cm^3/mol, ~, kg/m^3]
            
            % 2. Gas compressibility factor
            Z = P*V/(R*T_K);                                              % [-]
            
            % 3. Fugacity coefficients (Spycher et al., 2003)
            %   3.1 CO2
            fa = log(V/(V-b_co2));
            fb = b_co2/(V-b_co2);
            fc = (2*a_co2/(R*T_K^1.5*b_co2)) * log((V+b_co2)/V);
            fd = a_co2*b_co2/(R*T_K^1.5*b_co2^2) * (log((V+b_co2)/V) - b_co2/(V+b_co2));
            fe = log(Z);
            phi_co2 = exp(fa + fb - fc + fd - fe);                        % [-]
            
            %   3.2 Brine (with CO2)
            fbb = b_h2o/(V-b_co2);
            fcb = (2*a_h2o_co2/(R*T_K^1.5*b_co2)) * log((V+b_co2)/V);
            fdb = a_co2*b_h2o/(R*T_K^1.5*b_co2^2) * ...
                (log((V+b_co2)/V) - b_co2/(V+b_co2));                     % typo in Hassanzadeh et al., 2009
            phi_h2o = exp(fa + fbb - fcb + fdb - fe);                     % [-]
            
            % 4. CO2 molality in pure H2O at P, T conditions (m0_co2)
            B = phi_co2*P/(m_h2o*kap0_co2_g) * exp(-(P-p0)*Vp_co2/(R*T_K));
            A = kap0_h2o/(phi_h2o*P) * exp((P-p0)*Vp_h2o/(R*T_K));
            % y_i = mole fraction of component i in the gaseous phase
            y_h2o = (1-B)/(1/A-B);                                        % [-]
            % x_i = mole fraction of component i in the aqueous phase
            x_co2 = B*(1-y_h2o);                                          % [-]
            % Check values within range
            if x_co2 < 0, warning(['x_co2 = ' num2str(x_co2) '. Set to 0'])
                x_co2 = 0;
            elseif x_co2 > 1, warning(['x_co2 = ' num2str(x_co2) '. Set to 1'])
                x_co2 = 1;
            end
            if y_h2o < 0, warning(['y_h2o = ' num2str(y_h2o) '. Set to 0'])
                y_h2o = 0;
            elseif y_h2o > 1, warning(['y_h2o = ' num2str(y_h2o) '. Set to 1'])
                y_h2o = 1;
            end
            x_h2o = 1 - x_co2;                                            % [-]
            % molality
            m0_co2 = m_h2o*x_co2/x_h2o;                                   % [m]
            
            % 5. CO2 molality in saline solution at P, T conditions (m_co2)
            %   5.1 CO2 activity coefficient (Duan and Sun, 2003; as in Spycher
            %   & Pruess, 2005; typos in Hassanzadeh et al., 2009)
            gamma_co2 = activityCO2DS2003(T_K, P, m_io);
            
            %   5.2 CO2 molality
            m_co2 = m0_co2/gamma_co2;
            
            % 6. H2O and CO2 mole fractions and equilibrium ratios in aqueous
            % saline solution with dissolved CO2.
            X_co2  = m_co2/(m_co2 + m_h2o + sum(nu.*m_salt));             % [-] mole fraction
            if saltVar == 1
                X_salt = sum(m_salt)/(m_co2 + m_h2o + sum(m_salt));       % not ionized
            else
                X_salt = sum(x_salt);
            end
            X_salt_fi = sum(nu.*m_salt)/(m_co2 + ...
                                         m_h2o + sum(nu.*m_salt));        % [-] fully ionized
            X_h2o     = 1 - X_co2 - X_salt;                               % [-]
            X_solv    = 1 - X_co2;                                        % [-]
            Y_h2o     = A*(X_solv - X_salt_fi);                           % [-]
            Y_co2     = 1-Y_h2o;                                          % [-]
           %K_co2     = Y_co2/X_co2;                                      % [-] equilibrium rat.
           %K_h2o     = Y_h2o/X_h2o;                                      % [-]
           
            if iterate == 1
                res = abs(Y_h2o - Y0_h2o);
                Y0_h2o = Y_h2o;   
                Y0_co2 = Y_co2;
                if  res < 1e-6
                    break
                elseif it == maxits
                    error('Iterative loop did not converge')
                end
                it = it + 1;
            else
                break
            end 
        end
        disp(['Iterations at p=' num2str(P) 'bar: ' num2str(it)]);
               
        % 7. Compute Rs and Rv
        Rs = rhox_brine_s*X_co2/(rhox_co2_s*(1-X_co2));                   % [V CO2 Sm^3/ V br Sm^3] at sc.
        Rv = rhox_co2_s*Y0_h2o/(rhox_h2o_s*Y0_co2);                       % [V H2O Sm^3/ V CO2 Sm^3]
    else
        [V, ~, rho_co2] = rho_co2_fcn(T_K, P, a_m, b_m);                  % [cm^3/mol, ~, kg/m^3]
        Z = NaN;
        phi_co2 = NaN;
    end
    
    % 8. Molar masses and mass fractions
    Mw_aq   = Mw_salts*X_salt + Mw_h2o*X_h2o + Mw_co2*X_co2;              % [kg/mol] avg Molar mass
    Mw_solv = (Mw_salts*X_salt + Mw_h2o*X_h2o)/X_solv;
    Mw_gas  = Mw_co2*Y0_co2 + Mw_h2o*Y0_h2o;
    if saltVar == 1
        w_salt    = X_salt*(Mw_salts/Mw_aq);
    end
    w_co2 = X_co2*(Mw_co2/Mw_aq);                                         % [-] mass frac.
    v_h2o = Y0_h2o*(Mw_h2o/Mw_gas);
       
    % 9. Aqueous phase (CO2 saturated; with dissolved salts) density
    %    (Garcia, 2001) and Gaseous phase density
    rho_brine = rho_b_fcn(T_K, P, w_salt);                                % [kg/m^3, 1/kPa]
    rho_aq    = rho_aq_fcn(Mw_co2, Mw_solv, X_co2, X_solv, V_phi, rho_brine);
    rho_gas   = rho_gas_fcn(Mw_gas, V/10^6);
    
    % 10. FVF aqueous saline solution with CO2 (= B_aq) and gaseous
    % (CO2-rich) solution with H2O
    B_b   = rho_brine_s/rho_brine; 
    B_aq  = rho_brine_s/(rho_aq*(1 - w_co2));                             % [m^3/Sm^3]
    B_gas = rho_co2_s/(rho_gas*(1 - v_h2o));                              % [ " ]
    
    % 11. Viscosity of aqueous and gaseous solutions
    w_co2  = X_co2*(Mw_co2/Mw_aq);
    mu_aq  = mu_aq_fcn(T_K, P, sum(m_salt), w_co2)*10^3;                  % [cP]
    
    mu     = [mu_aq_fcn(T_K, P, 0, 0), mu_co2_fcn(T_K, rho_co2)];
    mu_gas = mu_gas_fcn([Y0_h2o, Y0_co2], [Mw_h2o, Mw_co2], mu)*10^3;     % [ " ]
    
    % Assign variables to fields in result table
    [t.gas.Z_co2(n), t.gas.phi_co2(n), t.gas.gamma_co2(n), ...
     t.gas.y_h2o(n), t.gas.Y_h2o(n), t.gas.rho_g_kgm3(n), ...
     t.gas.Rv_Sm3Sm3(n), t.gas.B_g_m3Sm3(n), t.gas.mu_g_cP(n)]  ...
     = deal(Z, phi_co2, gamma_co2, y_h2o, Y_h2o, rho_gas, Rv, B_gas, mu_gas);
           
    [t.aq.m_co2(n), t.aq.x_co2(n), t.aq.X_co2(n), t.aq.X_salt(n), ...
     t.aq.rho_b_kgm3(n), t.aq.rho_aq_kgm3(n), t.aq.Rs_Sm3Sm3(n), ...
     t.aq.B_b_m3Sm3(n), t.aq.B_aq_m3Sm3(n), t.aq.mu_aq_cP(n)]  ...
     = deal(m_co2, x_co2, X_co2, X_salt, rho_brine, rho_aq, Rs, B_b, B_aq, mu_aq);
end


%% Plots
if figs == 1
    col_co2  = [180, 0, 0; ...
                100, 100, 100; ...
                0, 99, 67]./255;
    col_brine  = [0, 71, 148]./255;
    %axisarg    = {'FontSize', 11, 'TickLabelInterpreter','latex'};
    latx       = {'Interpreter','latex'};
    fontsz_tit = {'fontsize',14};
    
    h1 = figure(1); % Compare with Fig. 8 in Spycher et al. (2003)
    yyaxis left
    plot(p,t.gas.Z_co2, 'color', col_co2(1,:)); grid on; ylim([0 1]);
    ylabel('$Z$ [-]','fontsize', 12, latx{:})
    yyaxis right
    plot(p,t.gas.phi_co2, 'color', col_co2(2,:), 'LineStyle', '--'); grid on;
    ylim([0 1]);
    ylabel('$\phi$ [-]','fontsize', 12, latx{:})
    xlabel('$p$ [bar]','fontsize', 12, latx{:})
    ax = gca;
    ax.YAxis(1).TickValues = 0:.1:1;
    ax.YAxis(2).TickValues = 0:.1:1;
    ax.YAxis(1).Color = col_co2(1,:); ax.YAxis(1).FontSize = 10;
    ax.YAxis(2).Color = col_co2(2,:); ax.YAxis(2).FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.XTick = [0 100 200 300 400 500 600];
    xlim([0 pmax]);
    legend('$Z$', '$\phi$', latx{:}, 'location', 'southwest', 'box', 'off')
    title(['CO$_2$ coefs., $T$=' num2str(T_C) '$^\circ$C'], latx{:}, ...
        fontsz_tit{:})
    set(h1, 'Position', [600, 600, 300, 260])
    
    h2 = figure(2); % Compare with Fig. 2 in Spycher and Pruess (2005)
    yyaxis left
    plot(p,t.aq.m_co2,'color', col_co2(1,:)); grid on; ylim([0 2]);
    ylabel('CO$_{2(aq)}$ [m]','fontsize', 12, latx{:})
    yyaxis right
    plot(p,t.gas.Y_h2o*1000,'color', col_brine, 'LineStyle', '--'); grid on; ylim([4 14]);
    ylabel('$\Psi_{\mathrm{H}_2\mathrm{O}}\times10^3$ [-]','fontsize', 12, ...
        latx{:})
    xlabel('$p$ [bar]','fontsize', 12, latx{:})
    ax = gca;
    ax.YAxis(1).TickValues = [0, 0.5, 1.0, 1.5, 2.0];
    ax.YAxis(2).TickValues = 4:14;
    ax.YAxis(1).Color = col_co2(1,:); ax.YAxis(1).FontSize = 10;
    ax.YAxis(2).Color = col_brine; ax.YAxis(2).FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.XTick = [0 100 200 300 400 500 600];
    xlim([0 pmax]);
    legend('CO$_{2\mathrm{(aq)}}$', '$\Psi_{\mathrm{H}_2\mathrm{O}}$', ...
        latx{:}, 'location', 'northeast', 'box', 'off')
    title(['CO$_{2\mathrm{(aq)}}$ and H$_2$O$_{\mathrm{(g)}}$, $T$=' ...
        num2str(T_C) '$^\circ$C'], latx{:},  fontsz_tit{:})
    set(h2, 'Position', [600, 600, 300, 260])
    
    h3 = figure(3); % Compare with Fig. 4 in Hassanzadeh et al. (2008)
    yyaxis left
    plot(p/10,t.aq.x_co2,'color', col_co2(1,:)); grid on; ylim([0.005 0.03]);
    ylabel('$\chi_{\mathrm{CO}_2}$ [-]','fontsize', 12, latx{:})
    yyaxis right
    plot(p/10,t.gas.y_h2o,'color', col_brine, 'LineStyle', '--'); grid on; ylim([0 0.03]);
    ylabel('$\psi_{\mathrm{H}_2\mathrm{O}}\times10^3$ [-]','fontsize', 12, ...
        latx{:})
    xlabel('$p$ [MPa]','fontsize', 12, latx{:})
    ax = gca;
    ax.YAxis(1).TickValues = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03];
    ax.YAxis(2).TickValues = [0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03];
    ax.YAxis(1).Color = col_co2(1,:); ax.YAxis(1).FontSize = 10;
    ax.YAxis(2).Color = col_brine; ax.YAxis(2).FontSize = 10;
    ax.XAxis.FontSize = 10;
    ax.XTick = [0 10 20 30 40 50 60 70 80];
    xlim([0 pmax/10]);
    legend('$\chi_{\mathrm{CO}_2}$', '$\psi_{\mathrm{H}_2\mathrm{O}}$', ...
        latx{:}, 'location', 'northeast', 'box', 'off')
    title(['CO$_{2\mathrm{(aq)}}$ and H$_2$O$_{\mathrm{(g)}}$, $T$=' ...
        num2str(T_C) '$^\circ$C'], latx{:},  fontsz_tit{:})
    set(h3, 'Position', [600, 600, 300, 260])
    
    h4 = figure(4); % Compare with Fig. 5 in Hassanzadeh et al. (2008). Typo in
    % the caption? it seems that T should be ~46C
    yyaxis left
    plot(p/10, t.aq.Rs_Sm3Sm3, 'color', col_co2(1,:));
    %ylim([0 20])
    ylabel('R$_s$ [Sm$^3$ CO$_2$/Sm$^3$ brine]','fontsize', 14, latx{:})
    %set(gca,'Ytick',[0 5 10 15 20])
    ax = gca;
    ax.YAxis(1).Color = col_co2(1,:); ax.YAxis(1).FontSize = 11;
    ytl = get(gca, 'YTick');
    yyaxis right
    plot(p/10, t.aq.B_aq_m3Sm3, 'color', col_brine, 'LineStyle', '--'); grid on;
    %ylim([1.012 1.032]); 
    %set(gca,'Ytick',[1.012 1.016 1.020 1.024 1.028 1.032])
    ylabel('B$_b$ [m$^3$/Sm$^3$]','fontsize', 14, latx{:})
    ax = gca;
    ax.YAxis(2).Color = col_brine; ax.YAxis(2).FontSize = 11;
    ytr = get(gca, 'YTick');
    ytrv = linspace(min(ytr), max(ytr), numel(ytl));
    ytrc = compose('%.3f', ytrv);
    set(gca, 'YTick', ytrv, 'YTickLabel', ytrc)
    
    xlabel('$p$ [MPa]','fontsize', 14, latx{:})
    xlim([0 50]);
    ax.XAxis.FontSize = 11;
    ax.XTick = [0 10 20 30 40 50];
    grid on; 
    legend('R$_s$', 'B$_b$', latx{:}, 'location', 'southeast', ...
           'box', 'off', 'fontsize', 12)
    title(['$T$=' num2str(T_C) 'C, $s$=' num2str(S{end}) S{1}], ...
        latx{:},  fontsz_tit{:})
    set(h4, 'Position', [600, 600, 400, 300])
    
    f5 = figure(5); % compare with Fig. 5b, d in Hassanzadeh et al. (2008)
    % subplot 1, density
    subplot(1,2,1)
    hold on
    p1 = plot(t.P_bar, t.aq.rho_b_kgm3, '--', 'color', col_brine, ...
        'linewidth', 1, 'DisplayName', '$\rho_\mathrm{b}$');
    p2 = plot(t.P_bar, t.aq.rho_aq_kgm3, '-', 'color', col_brine, ...
        'linewidth', 1, 'DisplayName', '$\rho_\mathrm{b}$, CO$_2$ sat.');
    p3 = plot(t.P_bar, t.gas.rho_g_kgm3, '-', 'color', col_co2(1,:), ...
        'linewidth', 1, 'DisplayName', '$\rho_\mathrm{g}$');
    hold off
    grid on
    xlabel('$p$ [bar]', latx{:}, 'fontsize', 12)
    ylabel('$\rho_\alpha$ [kg/m$^3$]', latx{:}, 'fontsize', 12)
    xlim([0 t.P_bar(end-2)])
    ylim([0 1200])
    legend([p1 p2, p3], latx{:}, 'fontsize', 10, 'location', 'southeast')
    % subplot 2, viscosity
    subplot(1,2,2)
    hold on
    p1 = plot(t.P_bar, t.aq.mu_aq_cP, '-', 'color', col_brine, ...
        'linewidth', 1, 'DisplayName', '$\mu_\mathrm{b}$, CO$_2$ sat.');
    p2 = plot(t.P_bar, t.gas.mu_g_cP, '-', 'color', col_co2(1,:), ...
        'linewidth', 1, 'DisplayName', '$\mu_\mathrm{g}$');
    hold off
    grid on
    xlabel('$p$ [bar]', latx{:}, 'fontsize', 12)
    ylabel('$\mu_\alpha$ [cP]', latx{:}, 'fontsize', 12)
    legend([p1 p2], latx{:}, 'fontsize', 10, 'location', 'southeast')
    set(gca,'yscale','log')
    ylim([1e-2 2])
    xlim([0 t.P_bar(end-2)])
    set(f5, 'Position', [200, 200, 600, 300])
end


%% Save PVT table
% Write to file
if nargin > 6  
    fileName_aq  = ['PVTO_T' num2str(T_K) '_P' num2str(p(1)) 'P' ...
                    num2str(pmax) '_S' s_type  num2str(s_val) s_unit];
    if max(t.gas.Rv_Sm3Sm3) > 0
        fileName_gas = ['PVTG_T' num2str(T_K) '_P' num2str(p(1)) 'P' ...
                        num2str(pmax) '_S' s_type  num2str(s_val) s_unit];
        tgas = table(t.P_bar, t.gas.Rv_Sm3Sm3, t.gas.B_g_m3Sm3, t.gas.mu_g_cP);
        tgas.Properties.VariableNames = {'P_bar', 'Rv_Sm3Sm3', 'B_g_m3Sm3', 'mu_g_cP'};
        disp('-----------------IMPORTANT NOTE----------------------------')
        disp('The last 2 rows in PVTG must be substituted by a single row.')
        disp('This single row should be for last (max) p value used as input and should be:')
        disp('Rv = 0, and corresponding Bgas and viscosity.')
        disp('To obtain, run model for iterate = 0 (PVDG) and get corresponding row.')
        disp('-----------------------------------------------------------')
    else
        fileName_gas = ['PVDG_T' num2str(T_K) '_P' num2str(p(1)) 'P' ...
                        num2str(pmax) '_S' s_type  num2str(s_val) s_unit];
        tgas = table(t.P_bar, t.gas.B_g_m3Sm3, t.gas.mu_g_cP);
        tgas.Properties.VariableNames = {'P_bar', 'B_g_m3Sm3', 'mu_g_cP'};
    end
    taq = table(t.aq.Rs_Sm3Sm3, t.P_bar, t.aq.B_aq_m3Sm3, t.aq.mu_aq_cP);
    taq.Properties.VariableNames = {'Rs_Sm3Sm3' 'P_bar' 'B_aq_m3Sm3' 'mu_b_cP'};
    writetable(taq, [directory fileName_aq '.txt'], 'Delimiter', 'tab');
    writetable(tgas, [directory fileName_gas '.txt'], 'Delimiter', 'tab');
end

end
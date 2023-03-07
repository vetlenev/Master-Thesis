%% Along-fault permeability based on single or multiple smears
%
% 
%
%
%
clear, close all
%% Preliminaries
% discontinuous smears model
ispar = false;                                   % true for alpha = 0

% Independent variables
clay_thickness  = 25;                           % m
fault_dip       = 30:5:90;                      % degrees
SSF             = 1:8;                          % shale smear factor
fault_throw     = SSF*clay_thickness;           % m
fTR             = [0.005 0.01 0.05 0.1];        % fault_thick/fault_disp
SSFc            = 4;                            % critical SSF

% Preallocate 
Ls      = cell(numel(fTR), 1);
k_dip   = cell(numel(fTR), 1);
Lf      = cell(numel(fTR), 1);
Psmear  = cell(numel(fTR), 1);
Td      = cell(numel(fTR), 1);
% Flag: 1: continuous (SSF <= SSFc); 
%       2: discontinuous (SSF > SSFc); 
%       3: continuous (due to smear thickness; SSF > SSFc);
flag    = cell(numel(fTR), 1);                  

% Compute
for j = 1:numel(fTR)
    for k = 1:numel(fault_throw)
        %% 1. Define material properties for sand and clay and geometry
        % Sand-based FZ material
        sand.phi    = 35;
        sand.theta  = 45 + sand.phi/2;
        sand.kzz    = 50;                        % along fault permeability [mD]
        
        % Clay
        clay.phi   = 5;
        clay.theta = 45 + clay.phi/2;
        clay.T     = 25;
        clay.dip   = 0;
        
        % Fault
        fault.t     = fault_throw(k);
        fault.dip   = fault_dip;
        fault.gamma = 90 - fault.dip;
        fault.D     = fault.t./sind(fault.dip);
        fault.T     = fTR(j)*fault.D;             % fault core
        
        % Fault-length of interest (includes source in both FW and HW)
        Lf{j}(k,:) = (fault.t + clay.T) ./ sind(fault.dip);
        
        %% 2. Smear in the fault
        % Permeability
        smear.kzz    = 1e-3;                      % along fault permeability [mD]
        % Angles
        smear.delta = atand( fault.t ./ (fault.T./sind(fault.dip) + ...
                             fault.t.*cotd(fault.dip)) );
        smear.alpha = 90 - fault.gamma - smear.delta;
        
        % Average smear thickness of Egholm et al. (2008). Here, the fault
        % dip is not equal to the sand theta (failure) angle (see Fig. 2 in the
        % paper); however, we make the assumption that synsedimentary
        % faults developed when at shallow depths and following MC's
        % failure theory, so that the amount of clay within the fault zone
        % is consistent with this. However, these faults can currently have
        % much gentler dips at depth, so that this is only strictly valid for
        % synsedimentary faults that do not transfer further displacement
        % at depth (where the fault no longer dips with the failure angle
        % of the sand).
        Tse = min([ (cotd(clay.theta) - cotd(sand.theta))*clay.T^2 ./ Lf{j}(k,:); ...
                     fault.T ]);
        
        % Average smear thickness
        %smear.Ts = min((clay.T/fault.t)*fault.T, Tse);  % minimum assumption
        smear.Ts = Tse;                                  % Egholm assumption (more realistic)
        
        
        % Average apparent clay thickness along the fault
        smear.Ls_avg = smear.Ts./sind(smear.alpha);
        
        % Apparent thickness of smear at the footwall cutoff
        smear.Tap = clay.T./(cosd(clay.dip).*sind(fault.dip));
        smear.LsA = smear.Tap;
        
        %% 3. Fault-along permeability
        if SSF(k) <= SSFc
            flag{j}(k, :) = ones(1,numel(fault_dip));
            Psmear{j}(k,:)= ones(1,numel(fault_dip));
            Ls{j}(k,:)    = min([smear.LsA; smear.Ls_avg]);
            ws            = Ls{j}(k,:) ./ Lf{j}(k,:);
            assert(all(ws <= 1));
            wss           = 1 - ws;
            k_dip{j}(k,:) = 1./(ws/smear.kzz + wss/sand.kzz);
       
        else % potentially discontinuous smears
            if ~ispar
                smear.Tap    = clay.T./(cosd(clay.dip).*sind(smear.delta));
                Ls{j}(k,:)   = smear.Tap .* SSFc;
                Psmear{j}(k,:) = Ls{j}(k,:) .* cosd(smear.alpha) ./ Lf{j}(k,:);
                smear.Td     = smear.Ts .* cosd(smear.alpha) + Ls{j}(k,:) .* sind(smear.alpha);
            else
                smear.Tap    = clay.T./(cosd(clay.dip).*cosd(fault.gamma));
                Ls{j}(k,:)   = smear.Tap .* SSFc;
                Psmear{j}(k,:) = Ls{j}(k,:) ./ Lf{j}(k,:);
                smear.Td     = smear.Ts;
            end
            
            % these are discontinuous
            disc             = smear.Td < fault.T;
            if ~ispar
                smear_free       = fault.T - Ls{j}(k,:) .* sind(smear.alpha);
                y                = smear.Ts/2 - smear_free./cosd(smear.alpha);
                disc2            = y; 
                disc2(disc2>0)   = false; % continuous due to smear thickness
                disc2(disc2<0)   = true;  % discontinuous
                disc             = any([disc; disc2] == true, 1);
            end
            flag{j}(k,disc)  = 2;     
            ws               = (smear.Td(disc) ./ fault.T(disc)) .* Psmear{j}(k,disc);
            assert(all(ws <= 1));
            wss              = 1 - ws;
            k_dip{j}(k,disc) = ws*smear.kzz + wss*sand.kzz;
            
            % continuous due to smear thickness
            cont                  = ~disc;
            if ~ispar
                L                 = y(cont) ./ sind(smear.alpha(cont));                       
                flag{j}(k,cont)   = 3;
                Ls{j}(k,cont)     = min([smear.LsA(cont); smear.Ls_avg(cont); L]);
                ws                = Ls{j}(k,cont) ./ Lf{j}(k,cont);
            else
                % consider length min + length min2 < Lf! when n > 1 
                ws                = Psmear{j}(k,cont);
            end            
            assert(all(ws <= 1));
            wss               = 1 - ws;
            k_dip{j}(k,cont)  = 1./(ws/smear.kzz + wss/sand.kzz);
        end
    end
end

%% Plots
colors = {[0.55 0.55 0.55],'k','c','r'};
dips   = [2 4 7 11];
latx   = {'interpreter', 'latex'};
axsz   = {'fontsize', 13};
lgsz   = {'fontsize', 11};

subplot(1,4,1)
semilogy([SSFc SSFc], [smear.kzz sand.kzz], '-', 'color', [0.8 0.8 0.8], 'linewidth', 2)
hold on
semilogy(SSF, k_dip{1}(:,dips(1)), 'o-','color', colors{1}, 'markersize', 9, 'markerFaceColor', colors{1})
semilogy(SSF, k_dip{1}(:,dips(2)), 's-','color', colors{2}, 'markersize', 8, 'markerFaceColor', colors{2})
semilogy(SSF, k_dip{1}(:,dips(3)), 'd-','color', colors{3}, 'markersize', 5, 'markerFaceColor', colors{3})
semilogy(SSF, k_dip{1}(:,dips(4)), '^-','color', colors{4}, 'markersize', 3, 'markerFaceColor', colors{4})
hold off
legend({'SSF$_\mathrm{c}$', ['$\beta =$ ' num2str(fault_dip(dips(1))) '$^\circ$'], ...
        [num2str(fault_dip(dips(2))) '$^\circ$'], [num2str(fault_dip(dips(3))) '$^\circ$'], ...
        [num2str(fault_dip(dips(4))) '$^\circ$']}, latx{:}, lgsz{:}, 'box', 'off', ...
        'location', 'northwest')
xlabel('SSF $= \frac{t}{T_\mathrm{c}}$', latx{:}, axsz{:})
xlim([0 SSF(end)])
xticks(SSF)
ylabel('$\hat{k_{zz}}$ [mD]', latx{:}, axsz{:})
ylim([smear.kzz sand.kzz])
yticks([10^-3 0.01 0.1 1 10 50])
ax = gca; 
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
title(['$\frac{\mathrm{f}_\mathrm{T}}{D} = $' num2str(fTR(1))], latx{:}, 'fontsize', 13)
grid on

subplot(1,4,2)
semilogy([SSFc SSFc], [smear.kzz sand.kzz], '-', 'color', [0.7 0.7 0.7], 'linewidth', 2)
hold on
semilogy(SSF, k_dip{2}(:,dips(1)), 'o-','color', colors{1}, 'markersize', 9, 'markerFaceColor', colors{1})
semilogy(SSF, k_dip{2}(:,dips(2)), 's-','color', colors{2}, 'markersize', 8, 'markerFaceColor', colors{2})
semilogy(SSF, k_dip{2}(:,dips(3)), 'd-','color', colors{3}, 'markersize', 5, 'markerFaceColor', colors{3})
semilogy(SSF, k_dip{2}(:,dips(4)), '^-','color', colors{4}, 'markersize', 3, 'markerFaceColor', colors{4})
hold off
xlim([0 SSF(end)])
xticks(SSF)
ylim([smear.kzz sand.kzz])
yticks([10^-3 0.01 0.1 1 10 50])
ax = gca; 
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
title(['$\frac{\mathrm{f}_\mathrm{T}}{D} = $' num2str(fTR(2))], latx{:}, 'fontsize', 13)
grid on

subplot(1,4,3)
semilogy([SSFc SSFc], [smear.kzz sand.kzz], '-', 'color', [0.7 0.7 0.7], 'linewidth', 2)
hold on
semilogy(SSF, k_dip{3}(:,dips(1)), 'o-','color', colors{1}, 'markersize', 9, 'markerFaceColor', colors{1})
semilogy(SSF, k_dip{3}(:,dips(2)), 's-','color', colors{2}, 'markersize', 8, 'markerFaceColor', colors{2})
semilogy(SSF, k_dip{3}(:,dips(3)), 'd-','color', colors{3}, 'markersize', 5, 'markerFaceColor', colors{3})
semilogy(SSF, k_dip{3}(:,dips(4)), '^-','color', colors{4}, 'markersize', 3, 'markerFaceColor', colors{4})
hold off
xlim([0 SSF(end)])
xticks(SSF)
ylim([smear.kzz sand.kzz])
yticks([10^-3 0.01 0.1 1 10 50])
ax = gca; 
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
title(['$\frac{\mathrm{f}_\mathrm{T}}{D} = $' num2str(fTR(3))], latx{:}, 'fontsize', 13)
grid on

subplot(1,4,4)
semilogy([SSFc SSFc], [smear.kzz sand.kzz], '-', 'color', [0.7 0.7 0.7], 'linewidth', 2)
hold on
semilogy(SSF, k_dip{4}(:,dips(1)), 'o-','color', colors{1}, 'markersize', 9, 'markerFaceColor', colors{1})
semilogy(SSF, k_dip{4}(:,dips(2)), 's-','color', colors{2}, 'markersize', 8, 'markerFaceColor', colors{2})
semilogy(SSF, k_dip{4}(:,dips(3)), 'd-','color', colors{3}, 'markersize', 5, 'markerFaceColor', colors{3})
semilogy(SSF, k_dip{4}(:,dips(4)), '^-','color', colors{4}, 'markersize', 3, 'markerFaceColor', colors{4})
hold off
xlim([0 SSF(end)])
xticks(SSF)
ylim([smear.kzz sand.kzz])
yticks([10^-3 0.01 0.1 1 10 50])
ax = gca; 
ax.XAxis.FontSize = 11;
ax.YAxis.FontSize = 11;
title(['$\frac{\mathrm{f}_\mathrm{T}}{D} = $' num2str(fTR(4))], latx{:}, 'fontsize', 13)
grid on


% % Compute Ls and Ws
% gamma = 90 - fault.dip;
% delta = atand( fault.t / (fault.T/sind(fault.dip) + ...
%                           fault.t*cotd(fault.dip)) );
% smear.alpha = 90 - gamma - delta;
% smear.LsB2  = smear.Tse/sind(smear.alpha); 
% 
% smear.LsB1 = ((clay.T/fault.t)*fault.T)/sind(smear.alpha);
% 
% % Correct for singularity as smear.alpha tends to 0;
% if smear.LsB2 > fault.D   
%     smear.alpha   = 0;
%     smear.LsB2    = fault.D;
% end
% if smear.LsB1 > fault.D
%     smear.LsB1 = fault.D;
% end
% smear.LsC  = clay.T/sind(fault.dip);
% 
% % Compute along-fault permeability
% % 1. Compute permeabilty along the bulk of the fault considering the 
% %    reduction in permeable space imposed by smear(s) located within 
% %    the fault (case A)
% ws        = smear.Tse/fault.T;
% fault.kyA = ws*clay.ky + (1-ws)*sand.ky;        % thickness-weighted arithmetic mean
% 
% % 2. Compute along-fault permeability, dominated by the portion of smear
% %    that the flow needs to cross (case C)
% smear.Ls = min([smear.LsB1, smear.LsB2, smear.LsC]);
% ws       = smear.Ls/fault.D;
% fault.ky = 1/(ws/clay.ky + (1-ws)/fault.kyA);  % length-weighted harm mean
%% Compute permeability of a fault with clay smears
%
% FW on the left, HW on the right
%
%
%
clc, clear, close all

% MRST load modules 
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui ...
           linearsolvers ls-proj ls-utils


%% 1. Initialize
% Input independent variables
% Note: fw and hw layers' thickness (T) is the true layer thickness
%       to agree with Egholm et a. (2008). Relation to the apparent thickness
%       is Tap = T/(cosBl*sinB) where Bl is layer dip angle and B is
%       fault dip angle.
f    = struct('T', 5, 'dip', 46.66);            % fThickness not from grid
sand = struct('perm', [100, 30]*milli*darcy, 'poro', 0.3, 'phi', 35);
clay = struct('perm', [0.01, 1e-4]*milli*darcy, 'poro', 0.05, 'phi', 10, 'SSFc', 3);
fw   = struct('T', repelem(25, 5), 'isclay', logical([1 0 1 0 1]), 'Bl', 0);
hw   = struct('T', [25 100], 'isclay', logical([1 0]), 'Bl', 0);
g    = struct('resL', 1, 'resT', 0.1);          % grid resolution [m]
inSmearL  = 'Lsmear';                           % Smear segment length [m] 
smears = 'objectLength';


%% 2. Create and populate grid
% Cartesian grid in MRST has indexing that starts at (0, 0) coordinate
% and runs faster on the X axis.
%
% Mapping grid and values for properties are stored in structure M.
fw.n = numel(fw.T); hw.n = numel(hw.T);
fw.id = 1:fw.n;
hw.id = fw.n:(fw.n+(hw.n-1));
assert(sum(hw.T) == sum(fw.T))

% 2.1 Create grid
f.t  = sum(hw.T(2:end));
f.D  = f.t./sind(f.dip);
Tap  = [fw.T ./ (cosd(fw.Bl)*sind(f.dip)) ...
        hw.T(2:end) ./ (cosd(hw.Bl)*sind(f.dip))];
f.L  = Tap + f.D;               % Length of interest
nels = max([round(f.T/g.resT), round(f.L(fw.n)/g.resL)]);
G = computeGeometry(cartGrid([nels,nels], [f.T f.L(fw.n)]));

% Plot grid
figure(1)
plotToolbar(G, G, 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1); 
axis equal; xlim([0 f.T]); ylim([0 f.L(fw.n)])
title(['Number of cells = ' num2str(G.cells.num)])


% 2.2 Properties

% 2.2.1 Compute Ts of each individual layer in FW and HW
clay.theta = 45 + clay.phi/2;
sand.theta = 45 + sand.phi/2;
f.Le = (f.t + hw.T(1)) / sind(sand.theta);      % Egholm et al. (2008) Lf
Ts = [(cotd(clay.theta) - cotd(sand.theta))*(fw.T.*fw.isclay).^2 ./ f.L(1:fw.n) ...
      (cotd(clay.theta) - cotd(sand.theta))*(hw.T.*hw.isclay).^2 ./ f.L(hw.id)];
Ts(fw.n+1) = [];

% 2.2.2 Diagonals that can have smear and those that cannot (only sand)
M.nDiagTot    = 2*G.cartDims(1) - 1;            % total number of diagonals
M.vals        = false(G.cartDims(1));           % actual matrix of 0s and 1s
M.units       = zeros(size(M.vals));            % Unit domain of each cell
M.unit        = 1:fw.n+hw.n-1;                  % Unit domain of each group
M.isclay      = [fw.isclay hw.isclay(2:end)];   % total units and clay or not  
M.nDiagMax    = round([G.cartDims(1)*(cumsum(fw.T)/sum(fw.T)) ...
                       (G.cartDims(1)-1)*(fliplr(cumsum(fliplr(hw.T)))/sum(hw.T))]); 
M.nDiagMax(fw.n) = sum(M.nDiagMax(fw.n:fw.n+1));
M.nDiagMax(hw.id(2)) = [];
        

% N Diagonals with potential smear
gamma = 90 - f.dip;
b     = f.T/sind(f.dip) + f.t*cotd(f.dip);
f.delta = atand(f.t/b);
f.alpha = 90 - gamma - f.delta;
f1    = atand(f.L ./ f.T);
zeta  = f1 - f.alpha;
Ls  = Ts ./ cosd(zeta);
fDiagL = sqrt(f.L .^2 + f.T^2); 
M.nDiag = min([round((Ls./fDiagL)*M.nDiagTot); M.nDiagMax]);
divLayerDiag = [round((M.nDiag(fw.n)-1)/2)+1 fix((M.nDiag(fw.n)-1)/2)];
if sum(M.nDiag) > M.nDiagTot
    M.nDiag(1:fw.n) = round(M.nDiag(1:fw.n) * ...
                      (G.cartDims(1)/sum([M.nDiag(1:fw.n-1) divLayerDiag(1)])));   
    M.nDiag(hw.id(2:end)) = round(M.nDiag(hw.id(2:end)) * ...
                            ((G.cartDims(1)-1)/sum([M.nDiag(hw.id(2:end)) divLayerDiag(2)])));
    divLayerDiag = [round((M.nDiag(fw.n)-1)/2)+1 fix((M.nDiag(fw.n)-1)/2)];
end

% N Diagonals with sand
clayFWDiag   = sum(M.nDiag(1:fw.n-1)) + divLayerDiag(1);
if  clayFWDiag > G.cartDims(1)
    vold       = clayFWDiag;
    clayFWDiag = G.cartDims(1);
    M.nDiag    = round(M.nDiag.*(clayFWDiag/vold)); 
end
sandFWDiag   = G.cartDims(1) - (clayFWDiag);
sandHWDiag   = max([(G.cartDims(1)-1)-(sum(M.nDiag)-clayFWDiag) 0]);
sandFWT      = fw.T(~fw.isclay); sandFWTtot = sum(sandFWT);
sandHWT      = hw.T(~hw.isclay); sandHWTtot = sum(sandHWT);
if ~M.isclay(fw.n)
    sandFWT(end) = sandFWT(end)/2; sandFWTtot = sum(sandFWT);
    sandHWT(1)   = sandHWT(1)/2;   sandHWTtot = sum(sandHWT);
end
M.nDiag(~fw.isclay) = round((sandFWT./sandFWTtot).*sandFWDiag);
M.nDiag(hw.id(~hw.isclay)) = M.nDiag(hw.id(~hw.isclay)) + ...
                             round((sandHWT./sandHWTtot).*sandHWDiag);
if sum(M.nDiag) < M.nDiagTot
   toadd = M.nDiagTot - sum(M.nDiag);
   disp(['Diagonals missing = ' num2str(toadd)])
   id = 1:toadd;
   M.nDiag(id) = M.nDiag(id) + 1;
elseif sum(M.nDiag) > M.nDiagTot
    tosub = [sum(M.nDiag(1:fw.n-1)) + divLayerDiag(1) - G.cartDims(1) ...
             divLayerDiag(2) + sum(M.nDiag(hw.id(2:end))) - (G.cartDims(1)-1)];
    disp(['Extra diagonals = ' num2str(tosub)])
    [~, idFW] = max(M.nDiag(1:fw.n));  
    [~, idHW] = max(M.nDiag(hw.id(2:end)));
    M.nDiag([idFW fw.n+idHW]) = M.nDiag([idFW fw.n+idHW]) - tosub;
end
assert(sum(M.nDiag) == M.nDiagTot) 

% Populate mapping matrix with all potential 1s and sure 0s
id1 = -(G.cartDims(1)-1);
idi = [id1 id1+cumsum(M.nDiag(1:end-1))];            % initial (all units)
idf = [idi(2:end)-1 (idi(end)-1)+M.nDiag(end)];      % final
for n = 1:numel(M.nDiag)
    if idf(n) < idi(n), idf(n) = idi(n); end    
    if M.isclay(n) == 1
        M.vals = full(spdiags(true(G.cartDims(1), abs(idi(n)-idf(n))+1), ...
                              -idf(n):-idi(n), M.vals));   
    else
        M.vals = full(spdiags(false(G.cartDims(1), abs(idi(n)-idf(n))+1), ...
                              -idf(n):-idi(n), M.vals));
    end   
    M.units = full(spdiags(n*ones(G.cartDims(1), abs(idi(n)-idf(n))+1), ...
                           -idf(n):-idi(n), M.units));
end
M.units = transpose(M.units); 
% M.vals is transposed at the end!

% 2.2.3 Smear continuity and location along possible diagonals (final 1s)
%       Here we assume that SSFc is the same in all layers.
% 2.2.3.1 Smear continuity
assert(abs(sum(Tap(1:fw.n)) - f.L(fw.n)) < 1e-5)
Tap2 = Tap;
Tap2(~M.isclay) = 0;
Ls      = Tap2*clay.SSFc;
Ds      = f.t/sind(f.delta);
M.Psmear  = Ls ./ Ds;                      % If Psmear < 1, discontinuous.
if any(M.Psmear > 1), M.Psmear(M.Psmear > 1) = 1; end

% 2.2.3.2 Smear distribution
if any(M.Psmear(M.isclay) < 1)
    idPsm = all([M.Psmear>0; M.Psmear<1]); 
    Psm = M.Psmear(idPsm);
    cunits = M.unit(idPsm);
    Lsmear = Ls(M.unit(idPsm));
    switch smears  
        case 'random'        % only constraint to be matched is P(smear)
            for j = 1:numel(cunits)
                ids = M.units==cunits(j);
                R = rand(size(M.vals(ids)));
                R(R>Psm(j)) = false;
                R(R>0) = true;
                M.vals(ids) = R;
            end
            
        case 'objectLength'         % place objects (segments) of length smear, and match P(smear)
            tolP = 1e-3;            % tolerance in deviation from P(smear)
            DiagCells = G.cartDims(1);
            cellSides    = unique(round(G.faces.areas, 4));
            assert(numel(cellSides) == 2);
            Dcell        = round(sqrt(sum(cellSides.^2)), 3);
            for j = 1:numel(cunits)
                disp(['Smear from source unit ' num2str(cunits(j)) ' is discontinuous.'])
                s = cunits(j);
                if ischar(inSmearL) && strcmp(inSmearL, 'Lsmear')
                    smearL = fix(Lsmear(j));
                elseif ischar(inSmearL)
                    error('If input as string, only Lsmear is accepted.')
                else
                    smearL = inSmearL;
                end
                claySegCellNum = round(smearL*(DiagCells/Ds));
                [subv, addv] = deal(round((claySegCellNum-1)/2), fix((claySegCellNum-1)/2));
                toadd        = [-subv:-1 0 1:addv];
                Dids         = [round((M.nDiag(cunits(j))-1)/2) ...
                                fix((M.nDiag(cunits(j))-1)/2)];
                Dvals_Mtest  = [round(G.cartDims(1)-Dids(1)):G.cartDims(1) ...
                                fliplr(fix(G.cartDims(1)-Dids(2)):G.cartDims(1)-1)];
                assert(numel(Dvals_Mtest) == sum(Dids)+1)
                % iterate to end up with number of cells with smear =
                % P(smear)
                smearCells = [];
                Pfin       = 0;
                maxIts     = 5;
                itNum      = 0;
                vdiv       = 1;
                sandCells  = 1:DiagCells;
                while itNum < maxIts
                    itNum          = itNum + 1;
                    Pfin0          = Pfin;                  
                    smearCellsT  = ((Psm(j)-Pfin)*DiagCells); 
                    smearNum     = round(smearCellsT/claySegCellNum);
                    if smearCellsT < 1
                        smearCellsT = round(smearCellsT);
                        if smearCellsT == 0
                            disp(['Tolerance of ' num2str(tolP) ...
                                  ' cannot be met: not enough cell resolution.'])
                        break
                        else
                            smearCellsT = 1;
                        end
                    end
                    if smearNum == 0
                        smearNum = 1;
                        claySegCellNum = round(smearCellsT);
                        [subv, addv] = deal(round((claySegCellNum-1)/2), fix((claySegCellNum-1)/2));
                        toadd        = [-subv:-1 0 1:addv];
                    end
                    for k = 1:smearNum
                        sandLim = find(diff(sandCells) > 1);
                        sandBounds = [1 sandLim+1; ...
                                      sandLim numel(sandCells)];
                        sandCellNum = diff(sandBounds)+1;
                        %sandSegNum = numel(sandBounds) - 1;
                        [~, sandSegId] = max(sandCellNum);
                        sandSegCells = sandCells(sandBounds(1,sandSegId)):sandCells(sandBounds(2,sandSegId));
                        sandSegCellNum = numel(sandSegCells);
                        if sandSegCellNum > claySegCellNum
                            smearOffset = round((rand(1)*(sandSegCellNum-claySegCellNum))) - ...
                                          round((sandSegCellNum-claySegCellNum)/2);
                        else
                            smearOffset = 0;
                        end
                        center = sandSegCells(1) + round(sandSegCellNum/2);
                        smearCenter = center + smearOffset;
                        smearCenter(smearCenter==0) = 1;
                        cellsToAdd     = repmat(smearCenter',1,claySegCellNum) + toadd;
                        cellsToAdd(cellsToAdd<1) = 1;
                        cellsToAdd(cellsToAdd>G.cartDims(1)) = G.cartDims(1);
                        % Updates
                        smearCells = unique([unique(cellsToAdd(:)); smearCells]);
                        sandCells = 1:DiagCells;
                        sandCells(smearCells) = [];  
                    end
                    %smearCenterIds = round((rand(smearNum, 1)*numel(sandCells)));
                    %smearCenterIds(smearCenterIds==0) = 1;
                    %smearCenter    = sandCells(smearCenterIds);
                    %cellsToAdd     = repmat(smearCenter',1,segmentCells) + toadd;
                    %cellsToAdd(cellsToAdd<1) = 1;
                    %cellsToAdd(cellsToAdd>G.cartDims(1)) = G.cartDims(1);
                    %smearCells     = unique([unique(cellsToAdd(:)); smearCells]);
                    dvals          = false(DiagCells,1);
                    dvals(smearCells) = true;
                    Mtest = transpose(full(spdiags(repmat(dvals, 1, M.nDiag(cunits(j))), ...
                                         [-Dids(1):-1 0 1:Dids(2)], false(DiagCells))));
                    Pfin = sum(sum(Mtest)) / sum(Dvals_Mtest);
                    %Pfin  = numel(smearCells)/DiagCells;
                    if abs(Psm(j) - Pfin) < tolP 
                        break
                    elseif Pfin > Psm(j)
                        disp('2nd Iterative loop required to match Psmear in the 2D subdomain:')
                        itNum2  = 0;
                        maxIts2 = 25;
                        while itNum2 < maxIts2
                            itNum2 = itNum2 + 1;
                            diffs  = diff(smearCells);
                            smearCellsOut = smearCells(diffs>1);
                            if numel(smearCellsOut) == 0
                                idRemove = round(rand(1));
                                idRemove(idRemove==0) = numel(smearCells);
                                smearCells(idRemove) = [];
                            else
                                idRemove = round(rand(1)*numel(smearCellsOut));
                                idRemove(idRemove==0) = 1;
                                smearCells(smearCells==smearCellsOut(idRemove)) = [];
                            end
                            dvals2             = false(DiagCells,1);
                            dvals2(smearCells) = true;
                            Mtest = transpose(full(spdiags(repmat(dvals2, 1, M.nDiag(cunits(j))), ...
                                                 [-Dids(1):-1 0 1:Dids(2)], false(DiagCells))));
                            Pfin2 = sum(Mtest(Mtest==true)) / sum(Dvals_Mtest);
                            if abs(Pfin2 - Psm(j)) < tolP
                                disp(['Tolerance of ' num2str(tolP) ...
                                      ' met in ' num2str(itNum2) ' iterations.'])
                                break
                            elseif itNum2 == maxIts2
                                disp('Tolerance cannot be met: max number of iterations reached')
                            elseif Pfin2 < Psm(j)
                                disp(['Tolerance of ' num2str(tolP) ...
                                      ' cannot be met: not enough cell resolution.'])
                                disp(['Number of iterations = ' num2str(itNum2)])
                                break
                            end
                        end
                        dvals = dvals2;
                        Pfin  = Pfin2;
                        break
                    end
                            
                end
                disp(['PSSF-derived P = ' num2str(Psm(j)) '. '...
                      'Final smear ' num2str(j) ' P = ' ...
                      num2str(Pfin) ' (' num2str(itNum) ' main iterations)'])
                disp(['Input smear length [m] = ' num2str(inSmearL) ...
                      '. Number of final smear segments = ' num2str(sum(diff(dvals)>0))]);
                M.vals = full(spdiags(repmat(dvals, 1, abs(idi(s)-idf(s))+1), ...
                                          -idf(s):-idi(s), M.vals));
                disp('___________________________________________________')
                %if cunits(j)==fw.n
                %    assert(sum(M.vals(M.units==fw.n)) == sum(Mtest(Mtest==true)))
                %end
            end   
            
        case 'variogram'    % spatially correlated
            
    end
    
else
    disp('Continuous smears: P(smear) = 1')
end
% To compromise spdiags truncation rules with smears chopped parallel
% to f.T and not to f.L, the last step is to transpose.
M.vals = transpose(M.vals);

% 2.2.4 Map diagonals to Grid indexing
%       Grid indexing starts at bottom left, columns (x) move faster. Matrix 
%       starts counting at top left. and rows (y) move faster.

% figure(9); spy(sparse(M.vals))       % check that smears are correctly placed in M.
idGrid = reshape(transpose(flipud(M.vals)), G.cells.num, 1);

% 2.2.5 Assign permeability and porosity
rock.perm = repmat(fliplr(sand.perm), [G.cells.num, 1]);
rock.perm(idGrid,1) = clay.perm(2);
rock.perm(idGrid,2) = clay.perm(1);

% Figure
latx = {'Interpreter', 'latex'};
hh = figure(2);
% subplot(1,2,1)
% plotToolbar(G, reshape(transpose(flipud(M.units)), G.cells.num, 1), ...
%            'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0.1); 
% xlim([0 f.T]); ylim([0 f.L(fw.n)]); c = colorbar;  set(c,'YTick', 1:max(M.unit))
% xlabel('$x$ [m]', latx{:}); ylabel('$y$ [m]', latx{:})
% title('Unit domains'); set(hh, 'position', [400, 0, 350, 1000]) 
% subplot(1,2,2)
colormap(copper)
plotToolbar(G, log10(rock.perm(:,2)/(milli*darcy)), 'EdgeColor', [0.2 0.2 0.2], 'EdgeAlpha', 0); 
xlim([0 f.T]); ylim([0 f.L(fw.n)]); c = colorbar; c.Label.Interpreter = 'latex';
c.Label.String = '$\log_{10} k_\mathrm{yy}$ [mD]'; c.Label.FontSize = 12; 
xlabel('$x$ [m]', latx{:}); ylabel('$y$ [m]', latx{:})
title('Permeability'); set(hh, 'position', [400, 0, 150, 1000]) 
%axis equal tight

%% 3. Simulations 
% GRAVITY? If yes, DIRECTION IN 2D (y) and must be rotated!  

% 4. Compute permeability
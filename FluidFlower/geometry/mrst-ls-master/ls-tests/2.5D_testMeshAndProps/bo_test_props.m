%% Black Oil Test: fluid properties
% Based on: singlePhaseAD.m + black oil examples
%
% The purposes of this test model are:
% 1. To use an external triangular mesh as input to generate an extruded 3D
%    mesh.
% 2. To assign fluid properties (equal in all mesh regions) using both the
%    ECLIPSE deck format and pre-defined anonymous functions. The fluids
%    are dry gas and live oil.
% 3. To create a well completion (gas injection) and use pressure boundary
%    conditions at the top and sides.
% 4. To run the model using the black oil AD solver simulateScheduleAD
%
clc, clear, close all force
mrstModule add ad-props ad-blackoil deckformat ad-core mrst-gui

%% Set up model geometry
meshpath = 'mrst-2019a/myprojects/2.5Dtestmesh/nodes_coordinates.dat';
vertices = dlmread(meshpath);
x = vertices(:,1); y = vertices(:, 2).*-1;
% delaunay triangulation and extrude grid
thick = [50 linspace(50,300,15)];
G = makeLayeredGrid(triangleGrid([x, y]), thick);
% Permute axes
G.nodes.coords = G.nodes.coords(:,[3 1 2]);
% Compute grid geometry
G = computeGeometry(G);
% Gravity
gravity reset on
% Plot Grid
plotGrid(G, 'FaceColor', 'k', 'FaceColor', [255 255 204]./255, 'FaceAlpha', 0.5, 'EdgeColor', [0.5 0.5 0.5]);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); set(gca,'Ydir','reverse')
view([-125, 25]); axis equal tight; 

%% Define rock & fluid model
% Poro and perm
rock = makeRock(G, 100*milli*darcy, 0.3); % homogeneous perm
%rock.perm(G.layerSize*round(G.numLayers/2)+1:end)=10*milli*darcy;

% Read data based on SPE1 gas and oil properties. For the deck, we create
% an ECLIPSE deck file which contains PROPS to use. The remaining
% keywords are just to avoid errors in reading the deck, but are not used
% here.
deck = readEclipseDeck('C:\Users\lsalo\matlab\mrst-2019a\myprojects\2.5Dtestmesh\eclipse_data_files\ftest.DATA');
deck = convertDeckUnits(deck);
fluid = initDeckADIFluid(deck);

% Plot data
figure(2) % Rel perms, dry gas visc and B, and rock compressibility
col = lines(2);
sg = linspace(0,fluid.krPts.g(3),25)';
p_plot = linspace(1,800,51)'*barsa;
name = {'dry gas','oil'};
subplot(1,3,1)
plot(sg, fluid.krG(sg), 'linewidth', 1, 'Color',col(1,:))
hold on
plot(sg, fluid.krOG(1-sg), 'linewidth', 1, 'Color',col(2,:)); hold off
grid on
legend(name{1},name{2},'Location','NorthWest');
xlabel('S_g [-]'); xlim([0 1])
title('Rel perms')
subplot(1,3,2)
yyaxis left
plot(p_plot/barsa, fluid.muG(p_plot), 'linewidth', 1, 'Color',col(1,:))
yyaxis right
plot(p_plot/barsa, 1./fluid.bG(p_plot), 'linewidth', 1, 'Color',col(2,:))
set(gca,'YScale','log')
grid on
legend({'\mu [kg/m*s]','B_G[-]'},'Location','NorthWest');
xlabel('p [bar]');
title('dry gas PV')
subplot(1,3,3)
plot(p_plot/barsa, fluid.pvMultR(p_plot), 'linewidth', 1, 'Color',col(1,:))
grid on
xlabel('p [bar]');
title('Pore V Mult.')
%plot(p/barsa, pv_r(1).*exp(cr*(p-p_r)),'LineWidth',2);

figure(3) % Live oil
pargs = {'LineWidth', 1, 'MarkerSize',5,'MarkerFaceColor',[.5 .5 .5]};
% Extract data from the input deck
pvto = deck.PROPS.PVTO{1};
rsd  = pvto.key([1:end end]);
pbp  = pvto.data([pvto.pos(1:end-1)' end],1);
Bod  = pvto.data([pvto.pos(1:end-1)' end],2);
muOd = pvto.data([pvto.pos(1:end-1)' end],3);
% "For Bo, we plot both the tabulated data as well as interpolated data at
% various combinations of dissolved gas and reservoir pressures. For each
% data point, we must first determine whether the state is saturated or not
% by comparing the actual amount of dissolved gas against the maximum
% possible amount of solved gas for this pressure." (showSPEfluids.m)
[M, N]    = deal(11,51);
[RsMax,pMax] = deal(max(rsd), max(pbp));
[rs,p_plot]    = meshgrid(linspace(10,RsMax-10,M), linspace(0,pMax,N));
Rv        = reshape(fluid.rsSat(p_plot(:)), N, M);
isSat     = rs >= Rv;
rs(isSat) = Rv(isSat);
Bo        = 1./reshape(fluid.bO(p_plot(:), rs(:), isSat(:)),N,M);
muO       = reshape(fluid.muO(p_plot(:), rs(:), isSat(:)),N,M);
% Make plots
subplot(1,3,1)
plot(pbp/barsa,rsd,'-o',pargs{:}); xlabel('Pressure [bar]'); ylabel('R_s [-]')
title('Solution gas-oil ratio')
% Formation-volume factor
subplot(1,3,2), hold on, %set(gca,'FontSize',14)
for j=1:M
    i = isSat(:,j);
    plot(p_plot(i,j)/barsa, Bo(i,j),'b-',p_plot(~i,j)/barsa,Bo(~i,j),'-r');
end
plot(pbp/barsa,Bod,'-bo',pargs{:});
hold off, axis tight, xlabel('Pressure [bar]'); ylabel('B_o [-]')
title('Oil formation-volume factor');
% Viscosity
subplot(1,3,3), hold on, %set(gca,'FontSize',14)
for j=1:M
    i = isSat(:,j);
    plot(p_plot(i,j)/barsa,  convertTo(muO(i,j),  centi*poise),'b-',...
         p_plot(~i,j)/barsa, convertTo(muO(~i,j), centi*poise),'-r');
end
plot(pbp/barsa, convertTo(muOd, centi*poise),'-bo',pargs{:});
hold off, axis tight, xlabel('Pressure [bar]'); ylabel('\mu_o [-]')
title('Oil viscosity');

%% Wells
simTime = 10*year;

numwell = 1;                                  % 4 horizontal injectors
injmass = 10^9*(simTime/year);                % 1 Mt/year of gas
injrate = injmass/(fluid.rhoGS*simTime*numwell);
%wellLoc = [25,2300,500; 25,2360,500; 25,2420,500; 25,2480,500];
wellLoc = [25,2360,500];
dirw = 'z';
distx = pdist2(G.cells.centroids(:,1), wellLoc(:,1));
disty = pdist2(G.cells.centroids(:,2), wellLoc(:,2));
distz = pdist2(G.cells.centroids(:,3), wellLoc(:,3));
dist = sqrt(distx.^2 + disty.^2 + distz.^2);
[~, cellInx] = min(dist,[],1);
W = addWell([ ], G, rock, cellInx, 'Name', 'I1', 'Dir', dirw, ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1]); % order is always 'WOG', here 'OG'

% 10y injection = simTime
timesteps = [linspace(0.01,0.01,10) linspace(0.1, 0.1, 9) linspace(1,1,9)]*year;

%% Assign model (Two phase live oil - dry gas)
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false); % to display

%% Initialize state
% Need to initialize 'p', 's', 'rs' and 'rv'
% p: impose vertical equilibrium (fully saturated in oil)
g = norm(gravity);
rho_wr = 1020*kilogram/meter^3;
p_r = g*rho_wr*50;                         % p(z0) below 50 m of water col.
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil

s0  = repmat([1 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [1 0] 'OG'
rs0 = zeros(G.cells.num, 1);            % no dissolved gas at the beginning
rv0 = 0;                                % dry gas
state0 = struct('s',s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);
%region = getInitializationRegionsBlackOil(model, 1, 'datum_pressure', p_r, 'rs', 0);
%state0 = initStateBlackOilAD(model, region);

%% Plot well and initial pressure
clf
show = true(G.cells.num,1);
show(cellInx) = false;
plotCellData(G,p0/barsa, show,'EdgeColor','k', 'FaceAlpha', 1);
c = colorbar; %xlim([0 5025]); ylim([0 5000]); zlim([0 1000]);
c.Label.String = 'p [bar]';
plotWell(G,W, 'height',10);
set(gca,'Xdir','reverse')
%set(gca,'Ydir','reverse')
view([55,25]), camproj perspective 
axis equal tight;

%% Boundary conditions

% add 'bc' to schedule (i.e. remove from line 144 and create here)
f = boundaryFaces(G);        
fp = f(G.faces.centroids(f,3) < 990);
%fval = G.faces.centroids(fp,2)*rho*norm(gravity);  
[z_0, z_max] = deal(min(G.faces.centroids(fp,3)), max(G.faces.centroids(fp,3)));
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
p_f_val = reshape(deval(equil, G.faces.centroids(fp,3)), [], 1);  clear equil
bc = addBC([], fp, 'pressure', p_f_val, 'sat', [1, 0]);  

%% Simulate schedule
% Set up the schedule containing both the wells and the timesteps
schedule = simpleSchedule(timesteps, 'W', W, 'bc', []);
%schedule = simpleSchedule(timesteps, 'W', W, 'bc', bc);

[wellSols, states, report] = ...
   simulateScheduleAD(state0, model, schedule);

figure(1)
for n=1:numel(states)
    states{n}.dp = states{n}.pressure - states{1}.pressure;  
end
id_seal = G.cells.centroids(:,3) < 150; % i.e. top "seal" layer
plotToolbar(G, states); set(gca,'Xdir','reverse')
view(55, 25); camproj perspective; axis equal tight 
plotWell(G,W); c = colorbar; caxis([0.7 1.1])

plotCellData(G,states{end}.rs, 'faceAlpha', 0.7);
set(gca,'Xdir','reverse')
view(55, 25); camproj perspective; axis equal tight 
c = colorbar; %C = flipud(load('colormap.txt')); colormap(C/255);
%caxis([0 90]); cmocean('thermal');
%% Test to check whether extruded grid axes change works
% Based on: singlePhaseAD.m
%
% Single-phase compressible AD solver
% The purpose of the example is to give the first introduction to how one
% can use the automatic differentiation (AD) class in MRST to write a flow
% simulator for a compressible single-phase model. For simplicity, the
% reservoir is assumed to be a rectangular box with homogeneous properties
% and no-flow boundaries. Starting from a hydrostatic initial state, the
% reservoir is produced from a horizontal well that will create a zone of
% pressure draw-down. As more fluids are produced, the average pressure in
% the reservoir drops, causing a gradual decay in the production rate.

clc, clear, close all force

% permute
permute = 1;
%% Set up model geometry
meshpath = 'mrst-2017b/myprojects/gom/2Dtestmesh/nodes_coordinates.dat';
vertices = dlmread(meshpath);
x = vertices(:,1); y = vertices(:, 2).*-1;
% delaunay triangulation and extrude grid
thick = [50 linspace(50,300,15)];
G = makeLayeredGrid(triangleGrid([x, y]), thick);
% Permute axes
if permute == 1
    G.nodes.coords = G.nodes.coords(:,[3 1 2]);
    %G.nodes.coords = G.nodes.coords(:,[1 3 2]); NOT A FUNCTIONAL GRID
    %G.cells.centroids = G.cells.centroids(:,[3 1 2]);
    idz = 3;
else
    idz = 2;
end
% Compute grid geometry
G = computeGeometry(G);

plotGrid(G, 'FaceColor', 'k', 'FaceColor', [255 255 204]./255, 'FaceAlpha', 0.5, 'EdgeColor', [0.5 0.5 0.5]);
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); 
set(gca,'Xdir','reverse')
if permute == 1
    view(3);
else
    view([150,77]);
end
axis equal tight; 

%% Define rock model
rock = makeRock(G, 100*milli*darcy, 0.3);

cr   = 1e-6/barsa;
p_r  = 200*barsa;
pv_r = poreVolume(G, rock);
pv   = @(p) pv_r .* exp( cr * (p - p_r) );

p = linspace(100*barsa,300*barsa,50);
plot(p/barsa, pv_r(1).*exp(cr*(p-p_r)),'LineWidth',2);


%% Define model for compressible fluid
mu    = 5*centi*poise;
c     = 1e-3/barsa;
rho_r = 850*kilogram/meter^3;
rhoS  = 750*kilogram/meter^3;
rho   = @(p) rho_r .* exp( c * (p - p_r) );

plot(p/barsa,rho(p),'LineWidth',2);

%% Assume a single horizontal well
if permute == 1
    wellLoc = [25,2300,500; 25,2360,500; 25,2420,500; 25,2480,500];
    dirw = 'y';
    %wellLoc = [2300,25,500; 2360,25,500; 2420,25,500; 2480,25,500];
    %dirw = 'x';
else
    wellLoc = [2300,500,25; 2360,500,25; 2420,500,25; 2480,500,25];
    dirw = 'x';
end
distx = pdist2(G.cells.centroids(:,1), wellLoc(:,1));
disty = pdist2(G.cells.centroids(:,2), wellLoc(:,2));
distz = pdist2(G.cells.centroids(:,3), wellLoc(:,3));
dist = sqrt(distx.^2 + disty.^2 + distz.^2);
[~, cellInx] = min(dist,[],1);
W = addWell([ ], G, rock, cellInx, 'Name', 'P1', 'Dir', dirw);

%% Impose vertical equilibrium
if permute == 1
    gravity reset on
else
    gravity reset y on
end
g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,idz)));
equil  = ode23(@(z,p) g .* rho(p), [z_0, z_max], p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,idz)), [], 1);  clear equil

%% Plot well and initial pressure
clf
show = true(G.cells.num,1);
show(cellInx) = false;
plotCellData(G,p_init/barsa, show,'EdgeColor','k', 'FaceAlpha', 0.5);
c = colorbar; %xlim([0 5025]); ylim([0 5000]); zlim([0 1000]);
c.Label.String = 'p [bar]';
plotWell(G,W, 'height',10);
set(gca,'Xdir','reverse')
%set(gca,'Ydir','reverse')
view([55,25]), camproj perspective 
axis equal tight;
pause

%% Compute transmissibilities
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);                          % Interior neighbors
hT = computeTrans(G, rock);                 % Half-transmissibilities
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
T  = T(intInx);                             % Restricted to interior

%% Define discrete operators
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
grad = @(x)C*x;
div  = @(x)-C'*x;
avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
spy(C)

%% Define flow equations
gradz  = grad(G.cells.centroids(:,idz));
v      = @(p)  -(T/mu).*( grad(p) - g*avg(rho(p)).*gradz );
presEq = @(p,p0,dt) (1/dt)*(pv(p).*rho(p) - pv(p0).*rho(p0)) ...
                      + div( avg(rho(p)).*v(p) );

%% Define well equations
wc = W(1).cells; % connection grid cells
WI = W(1).WI;    % well-indices
dz = G.cells.centroids(cellInx,idz);
%dz = W(1).dZ;    % connection depth relative to bottom-hole ATTENTION THIS
% IS WRONG FOR THE PERMUTED CASE [3 1 2]

p_conn  = @(bhp)  bhp + g*dz.*rho(bhp); %connection pressures
q_conn  = @(p,bhp) WI .* (rho(p(wc)) / mu) .* (p_conn(bhp) - p(wc));

rateEq = @(p,bhp,qS)  qS-sum(q_conn(p, bhp))/rhoS;
ctrlEq = @(bhp)       bhp-100*barsa;

%% Initialize for solution loop
[p_ad, bhp_ad, qS_ad] = initVariablesADI(p_init, p_init(wc(1)), 0);
nc = G.cells.num;
[pIx, bhpIx, qSIx] = deal(1:nc, nc+1, nc+2);

numSteps = 36;                  % number of time-steps
totTime  = 3*year;              % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-6;                % Newton tolerance
maxits   = 10;                  % max number of Newton its

sol = repmat(struct('time',[],'pressure',[],'bhp',[],'qS',[]),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', double(p_ad), ...
   'bhp', double(bhp_ad), 'qS', double(qS_ad));

%% Main loop
t = 0; step = 0;
hwb = waitbar(t,'Simulation ..');
while t < totTime,
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));
   % Newton loop
   resNorm = 1e99;
   p0  = double(p_ad); % Previous step pressure
   nit = 0;
   while (resNorm > tol) && (nit <= maxits)
      % Add source terms to homogeneous pressure equation:
      eq1     = presEq(p_ad, p0, dt);
      eq1(wc) = eq1(wc) - q_conn(p_ad, bhp_ad);
      % Collect all equations
      eqs = {eq1, rateEq(p_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)};
      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      % Update variables
      p_ad.val   = p_ad.val   + upd(pIx);
      bhp_ad.val = bhp_ad.val + upd(bhpIx);
      qS_ad.val  = qS_ad.val  + upd(qSIx);

      resNorm = norm(res);
      nit     = nit + 1;
      fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
   end

   if nit > maxits,
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, 'pressure', double(p_ad), ...
                            'bhp', double(bhp_ad), 'qS', double(qS_ad));
      waitbar(t/totTime,hwb);
   end
end
close(hwb);
%% Plot production rate and pressure decay
clf,
[ha,hr,hp] = plotyy(...
   [sol(2:end).time]/day, -[sol(2:end).qS]*day, ...
   [sol(2:end).time]/day, mean([sol(2:end).pressure]/barsa), 'stairs', 'plot');
set(ha,'FontSize',16);
set(hr,'LineWidth', 2);
set(hp,'LineStyle','none','Marker','o','LineWidth', 1);
set(ha(2),'YLim',[50, 300],'YTick',50:50:300);
xlabel('time [days]');
ylabel(ha(1), 'rate [m^3/day]');
ylabel(ha(2), 'avg pressure [bar]');
%pause

%% Plot pressure evolution
clf;
steps = [2 12 24 36];
for i=1:4
   subplot(2,2,i);
   plotCellData(G, (sol(steps(i)).pressure - p_init)/barsa, show,'EdgeColor',.5*[1 1 1]);
   plotWell(G,W,'color','k');
   set(gca,'Xdir','reverse');
   view([55,25]), %camproj perspective
   %caxis([115 205]);
   axis equal tight off; %zoom(1.4)
   %xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]'); 
   text(3000,100,-300,[num2str(round(steps(i)*dt/day)) 'day'],'FontSize',14);
end
h=colorbar('horiz','Position',[.1 .05 .8 .025]);
colormap(jet(55));
pause

%% Plot subrange of cells
show = false(G.cells.num,1);
ids = linspace(0,G.numLayers-1,16); ids = transpose(ids*G.layerSize);
idcells = cellInx+ids; idcells = reshape(idcells,numel(idcells),1);
show(idcells) = true;
plotCellData(G, (sol(steps(i)).pressure - p_init)/barsa, show,'EdgeColor',.5*[1 1 1]);
plotWell(G,W,'color','k');
set(gca,'Xdir','reverse');
view([55,25])
axis equal tight
text(3000,100,-300,[num2str(round(steps(i)*dt/day)) 'day'],'FontSize',14);
%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

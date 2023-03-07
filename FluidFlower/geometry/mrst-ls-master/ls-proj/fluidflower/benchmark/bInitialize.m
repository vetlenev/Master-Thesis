function [state0, p_r] = bInitialize(G, fluid, plts, makePlot)
%
%
%
g = norm(gravity);
rho_wr = fluid.rhoOS*kilogram/meter^3;
water_column = 1.5 - max(G.cells.centroids(:,3)); % m
p_r = 1*barsa + g*rho_wr*water_column; % p at shallowest z
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* fluid.bO(p,0,false)*fluid.rhoOS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
rv0 = zeros(G.cells.num, 1);             % no vaporized water at the beginning
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

if makePlot == 1
    %plot
    plts.fig3D(); plotCellData(G, p0/barsa, 'edgealpha', 0.2)
    plts.setAxProps(gca), colormap(jet), c = colorbar;
    axis equal off
end


end
function state0 = init_bo(G, f, water_column, ucids, figs)
%
%
%
g = norm(gravity);
rho_wr = f.rhoOS*kilogram/meter^3;
if isfield(ucids, 'upper') && ucids.upper == 2  % model does not reach the seabed
    p_r = load('ls-proj/gom/adbo_2.5D/data_files/p0TopAmp.mat'); % from another run with mesh that reaches the seabed
    assert(abs(min(G.cells.centroids(:,3)) - p_r.z1) < 1e-6)
    p_r = p_r.p1;
else
    p_r = g*rho_wr*water_column;  % model reaches the seabed, and we add the specified water column                
end
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* f.bO(p,0,false)*f.rhoOS, [z_0, z_max], p_r);
p0 = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil

s0  = repmat([1, 0], [G.cells.num, 1]);  % s: fully saturated in oil --> [0 1 0] if 'WOG'; [1 0] if 'OG'
rs0 = zeros(G.cells.num, 1);             % no dissolved gas at the beginning
rv0 = 0;                                 % dry gas
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

if figs == 1
    h = figure(40);
    cmap = flipud(cmap_agu());
    colormap(cmap);
    plotToolbar(G, p0/barsa), set(gca,'Xdir','reverse');
    view(55, 25); camproj perspective; axis equal off
    title('Initial pressure')
    c = colorbar;
    c.Label.String = '$p_0$ [bar]';  c.Label.FontSize = 22;
    c.Label.Interpreter = 'latex';   c.FontSize = 18;
    set(h, 'Position', [500, -400, 5000, 1000])
    
    h = figure(41);
    cmap = flipud(cmap_agu());
    colormap(cmap);
    plotToolbar(G, p0/barsa, [ucids{1} ucids{58}]), set(gca,'Xdir','reverse');
    view(55, 25); camproj perspective; axis equal off
    title('Initial p, reservoir and seal')
    c = colorbar;
    c.Label.String = '$p_0$ [bar]';  c.Label.FontSize = 22;
    c.Label.Interpreter = 'latex';   c.FontSize = 18;
    set(h, 'Position', [500, -400, 5000, 1000])
    
    figure(42)
    plotToolbar(G, s0(:,1)), set(gca,'Xdir','reverse');
    view(55, 25); camproj perspective; axis equal off
    c = colorbar; caxis([0 1]); title('Brine saturation')
end

end
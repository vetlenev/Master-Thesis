function plotConnected(C, Vcl, dim)
%
%
%
N = size(C, 1);
f(1, :) = sum(C(:, :, 1), 1);
f(2, :) = sum(C(:, :, 2), 1);
if nargin > 2 && dim == 3
    f(3, :) = sum(C(:, :, 3), 1);
    y = true;
else
    y = false;
end

rr = [255, 125, 125]/255;
bb = [125, 125, 255]/255;
latx = {'Interpreter', 'latex'};
fsz_ax = {'fontsize', 12};

fh = figure(21);
hold on
ncomp = size(f,1);
plot(repelem(Vcl(1), ncomp), f(:, 1)./N, '-k')
plot(repelem(Vcl(2), ncomp), f(:, 2)./N, '-k')
plot(repelem(Vcl(3), ncomp), f(:, 3)./N, '-k')
plot(repelem(Vcl(4), ncomp), f(:, 4)./N, '-k')
px = plot(Vcl, f(1,:)./N, 'sk', 'markerfacecolor', [0.4 0.4 0.4], ...
          'markersize', 8, 'displayname', 'x');
pz = plot(Vcl, f(2, :)./N, 'ob', 'markerfacecolor', bb, ...
          'markersize', 8, 'displayname', 'z');
if y
    py = plot(Vcl, f(3, :)./N, 'dr', 'markerfacecolor', rr, ...
              'markersize', 8, 'displayname', 'y');
end
grid on
hold off
xlabel('$\overline{V_\mathrm{cl}}$ [-]', fsz_ax{:}, latx{:})
ylabel('$\nu/N$ [-]', fsz_ax{:}, latx{:})
if y
    legend([px, py, pz], 'fontsize', 10, latx{:})
else
    legend([px, pz], 'fontsize', 10, latx{:})
end
set(fh, 'position', [200, 200, 350, 275]);
xlim([0.1 0.6]), xticks(0.1:.1:.6)
ylim([0 1]), yticks(0:.2:1);

end
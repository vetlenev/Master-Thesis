
%
%
%
clear, close all

% Input params
Fvsh = log10(logspace(0.001,1,50));
ksh  = 0.005;
kss  = 50;
Kcat = [0 0.01 0.1 0.4 0.7 1];
Cf   = [0 0.1 0.25 0.35 0.5 0.65 0.75 0.9 1];

% Compute & plot
ecf  = (1-2*Cf)';
a    = Fvsh.*ksh.^ecf;
b    = cell(1,numel(Kcat));
kf   = cell(1,numel(Kcat));
cols = hsv(numel(Cf));
leg  = {['$C_\mathrm{f}$ = ' num2str(Cf(1))], ['$C_\mathrm{f}$ = ' num2str(Cf(2))], ...
        ['$C_\mathrm{f}$ = ' num2str(Cf(3))], ['$C_\mathrm{f}$ = ' num2str(Cf(4))], ...
        ['$C_\mathrm{f}$ = ' num2str(Cf(5))], ['$C_\mathrm{f}$ = ' num2str(Cf(6))], ...
        ['$C_\mathrm{f}$ = ' num2str(Cf(7))], ['$C_\mathrm{f}$ = ' num2str(Cf(8))], ...
        ['$C_\mathrm{f}$ = ' num2str(Cf(9))]};
latx = {'Interpreter', 'latex'};
sz   = {'fontSize', 11};
for n = 1:numel(Kcat)
      b{n} = (1-Fvsh).*(kss*Kcat(n)).^ecf;
      kf{n}  = (a + b{n}).^ecf;
    
    % Figure
    subplot(2,3,n)
    h = semilogy(Fvsh, kf{n}, '-', 'linewidth', 1); set(h, {'color'}, num2cell(cols,2));
    if n == 1 
        legend(leg, latx{:}); 
        xlabel('$F_\mathrm{vsh}$', latx{:}, sz{:})
        ylabel('$k_\mathrm{f}$ [mD]', latx{:}, sz{:}); 
    end
    ylim([10^-5, 10^2])
    grid on
    title(['$K_\mathrm{cat}$ = ' num2str(Kcat(n))], latx{:}, 'fontSize', 11)
end
subplot(2,3,1)
title(['K$_\mathrm{cat}$=' num2str(Kcat(1)), ', $k_\mathrm{sh}$=' num2str(ksh), ...
       ', $k_\mathrm{ss}$=' num2str(kss) ' [mD]'], latx{:}, 'fontSize', 11)

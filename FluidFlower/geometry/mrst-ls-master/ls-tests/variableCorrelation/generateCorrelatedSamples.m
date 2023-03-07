%% Generate correlated samples using copulas
% This is an example showing how to generate dependent (correlated) random
% samples using a gaussian copula and custom marginal distributions other 
% than the normal.
%
% Requirements:
% Statistics and Machine Learning Toolbox
%
% Reference:
% https://www.mathworks.com/help/stats/copulas-generate-correlated-samples.html
%

clc; close all; clear

% 1. Marginal distribution ranges and parameters
phi.range = [8 16];
phi.distro = 'beta';
phi.param = [3 5];   % a and b

ssfc.range = [3 7];
ssfc.distro = 'uniform';

% 2. Generate pairs of correlated values using a Gaussian copula
n = 1000;
p = -0.7;                                   % linear correlation coeff.
U = copularnd('Gaussian', [1 p; p 1], n);

figure(1)
scatterhist(U(:,1), U(:,2), 'Marker', '.', 'color', 'k', 'Direction', 'out')
title(['$\rho =$ ' num2str(p)], 'Interpreter', 'latex', 'fontsize', 14)
xlabel('U_1')
ylabel('U_2')
xlim([0 1]); ylim([0 1])

% 3. Transform bivariate data to desired marginals and ranges
X = [betainv(U(:,1), phi.param(1), phi.param(2)),  U(:, 2)];    % marginals
phi.vals  = phi.range(1) + diff(phi.range) .* X(:, 1);          % range X1
ssfc.vals = ssfc.range(1) + diff(ssfc.range) .* X(:, 2);        % range X2

figure(2)
scatterhist(phi.vals, ssfc.vals, 'Marker', '.', 'color', [0.5 0.5 0.5], ...
            'Direction', 'out')
grid on
xlabel('$\phi$ [deg]', 'Interpreter', 'latex', 'fontsize', 14)
ylabel('SSF$_\mathrm{c}$ [-]', 'Interpreter', 'latex', 'fontsize', 14)
xlim(phi.range); ylim(ssfc.range)


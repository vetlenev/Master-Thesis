function [f, sand, clay, fw, hw, smear] = smearOrganizeIn(f, sand, clay, fw, hw, smear)
%
% Organize inputs and fill in footwall, hangingwall, fault and thickness
% variables from inputs given by the user.
%

% fill in fw, hw and f structures
fw.n = numel(fw.T); hw.n = numel(hw.T);
fw.id = 1:fw.n;
hw.id = fw.n+1:fw.n+hw.n;
assert(sum(hw.T) == sum(fw.T))
f.t  = sum(hw.T);
f.D  = f.t./sind(f.dip);

% Compute Ts of each individual layer in FW and HW
numClays = sum(fw.isclay) + sum(hw.isclay);
if numel(clay.phi) < numClays
    clay.phi = repelem(clay.phi, numClays);
end
clay.theta = 45 + clay.phi/2;
sand.theta = 45 + sand.phi/2;
clay.Tap  = [fw.T(fw.isclay) ./ (cosd(fw.Bl)*sind(f.dip)) ...
             hw.T(hw.isclay) ./ (cosd(hw.Bl)*sind(f.dip))];
f.L  = clay.Tap + f.D;                          % Current Lf 
f.Le = (f.t + [fw.T hw.T]) ./ sind(sand.theta); % Egholm et al. (2008) Lf (after fault formation)
termCot = cotd(clay.theta) - cotd(sand.theta);
fwId = 1:sum(fw.isclay);
hwId = 1+sum(fw.isclay):1+sum(fw.isclay)+(sum(hw.isclay) - 1);
smear.Ts = [termCot(fwId).*(fw.T(fw.isclay)).^2 ./ f.L(fwId) ...
            termCot(hwId).*(hw.T(hw.isclay)).^2 ./ f.L(hwId)];
        
% sand Tap
sand.Tap  = [fw.T(~fw.isclay) ./ (cosd(fw.Bl)*sind(f.dip)) ...
             hw.T(~hw.isclay) ./ (cosd(hw.Bl)*sind(f.dip))];

% Compute Ls of each smear (length of smears in the direction parallel to a
% diagonal from the top left or top right of fault to bottom right or
% bottom left of fault.
gamma = 90 - f.dip;
b     = f.T/sind(f.dip) + f.t*cotd(f.dip);
f.delta = atand(f.t/b);                    % true values of actual shear zone
f.alpha = 90 - gamma - f.delta;            % "                      % "
f1      = atand(f.D ./ f.T);
zeta  = f1 - f.alpha;                           
smear.Ls  = smear.Ts ./ cosd(zeta);        % if z is very close to 90, this is > than true.

end
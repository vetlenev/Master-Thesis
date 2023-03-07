function [M, smear] = smearContinuity(fw, f, M, sand, clay, smear)
%
% Determine whether smears are continuous or not and update M and smear
% structures.
% For now we assume that SSFc is the same in all layers.
%
assert(abs(sum(clay.Tap(1:sum(fw.isclay))) + sum(sand.Tap(1:sum(~fw.isclay))) - f.D) < 1e-5)
smear.L  = (clay.Tap/cosd(f.alpha)).*clay.SSFc;
smear.Ds = f.t/sind(f.delta);
M.Psmear = smear.L ./ smear.Ds;                 % If Psmear < 1, discontinuous.
if any(M.Psmear > 1) 
    M.Psmear(M.Psmear > 1) = 1; 
end

%M.Psmear(M.Psmear==0) = [];
%smear.L(smear.L==0) = [];
M.isclay = logical(M.isclay);
if numel(M.Psmear) > sum(M.isclay)              % Some smears overprinted by wider ones.
    M.Psmear(M.idSmearInRemoved) = [];
    smear.L(M.idSmearInRemoved) = [];
end

end
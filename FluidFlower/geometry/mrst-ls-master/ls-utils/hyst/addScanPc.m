function f = addScanPc(f, regi)
%
% February 2020: Support for capillary pressure hysteresis in the
%               nonwetting phase using Killough's (1976) model. This uses
%               the Land (1968) model to determine the trapped gas
%               saturation (Sgt) and then Killough's model to determine the
%               oil-gas capillary pressure along a scanning curve.
%
%               It requires that the bounding imbibition curve already be
%               in the fluid structure as f.pcOGib. Currently, this
%               can be obtained from assignSGOF (keyword SGFN in ECLIPSE
%               input file). Assumes that gas is the nonwetting phase.
%
% Compute scanning capillary pressure curves.
% sg  = current gas saturation
% sgi = maximum gas saturation achieved in the current run.
% __________________________________________________________
regi  = unique(regi);
nregi = numel(regi);

% Maximum and minimum gas saturation
if isfield(f, 'sWcon')
    sgmax = 1 - f.sWcon;
    if numel(sgmax) ~= nregi  % temporary
        sgmax = ones(1, nregi)*sgmax;
    end
    sgmin = zeros(1, nregi);
    warning('Min gas saturation (primary drainage) set to 0')
elseif isfield(f, 'krPts')
    f.sgtmax   = f.krPts.g(regi,2)';
    sgmax      = f.krPts.g(regi, 3)';
    sgmin      = f.krPts.g(regi, 1)';
else
    sgmax = ones(1, nregi);
    sgmin = zeros(1, nregi);
    warning('Max gas saturation (primary drainage) set to 1.')
    warning('Min gas saturation (primary drainage) set to 0')
end
assert(numel(sgmax) == nregi);

% Compute scanning curve for Pc and add to fluid object
pci = cell(1, nregi);
for n = 1:nregi
    pciVars = {sgmax(n), sgmin(n), f.pcOG{n}, f.pcOG{n+nregi}, f.sgtmax(n)};
    pci{n}  = @(sg, sgi) scanPc(sg, sgi, pciVars, f.ehystr);
end
f.pcOGi   = pci;

% ----------------- Compute Scanning Pc curve --------------------
    function [pc, tol, minSat, sgt] = scanPc(s, si, pciVars, opts)  %LS
        %
        % February 2020: Support for capillary pressure hysteresis in the
        %               nonwetting phase using Killough's (1976) model. This uses
        %               the Land (1968) model to determine the trapped gas
        %               saturation (Sgt) and then Killough's model to determine the
        %               oil-gas capillary pressure along a scanning curve.
        %
        %               It requires that the bounding imbibition curve already be
        %               in the fluid structure as f.pcOGib. Currently, this
        %               can be obtained from assignSGOF (keyword SGFN in ECLIPSE
        %               input file). Assumes that gas is the nonwetting phase.
        %
        % Compute scanning capillary pressure curves.
        % sg  = current gas saturation
        % sgi = maximum gas saturation achieved in the current run.
        % __________________________________________________________
        nci = numelValue(s);
        
        % Set tolerance for flow reversal and initialize
        assert(nci == numelValue(si))   % same n for both sg and sgi
        tol  = 1e-3;                     % must be the same as for kr hysteresis!
        minSat = 0.05;
        pc  = nan(nci, 1);
        AD = getSampleAD(s, si);
        if isa(AD, 'GenericAD')
            pc = double2GenericAD(pc, AD);
        elseif isa(AD, 'ADI')
            pc = double2ADI(pc, AD);
        end
        
        % Inputs
        sgmx      = pciVars{1};
        sgmn      = pciVars{2};
        pcOG      = pciVars{3};
        pcOGib    = pciVars{4};
        sgtmx     = pciVars{5};
        
        % Index for cells where scanning curves are computed. We set a
        % minimum of 5% gas saturation to activate hysteresis.
        imb  = all([value(s) + tol < value(si), ...
            value(si) > minSat+sgmn], 2);  % MUST BE THE SAME USED IN KR HYST
        
        if any(imb)
            pc(~imb) = pcOG(s(~imb));
            sgiv  = si(imb);
            sgv   = s(imb);
            
            % Killough's hysteresis model
            a     = opts{4};        % modif. param. for improved converg.
            A     = 1 + a*(sgmx - sgiv);
            C     = 1/(sgtmx - sgmn) - 1/(sgmx-sgmn);
            sgt   = sgmn + (sgiv-sgmn)./(A + C*(sgiv-sgmn));
            e     = opts{1};
            if strcmp(opts{6}, 'RETR')
                F =  (1./(sgiv-sgv+e) - 1/e) ./ ...
                    (1./(sgiv-sgt+e) - 1/e);
                F(F>1) = 1;                               % no negative Pc
                v = pcOG(sgv) + F.*(pcOGib(sgv) - pcOG(sgv));
            else
                error(['Support for new scanning curves for second' ...
                    ' flow reversal onwards not coded yet'])
            end
            pc(imb) = v;
            
            isAbove = pc > pcOG(value(s));
            if any(isAbove)
                pc(isAbove) = pcOG(value(s(isAbove))) - 1e-5;
            end
            
            % Checks
            pciv = value(pc);
            assert(all(~isnan(pciv)));                   % no nan values
            assert(all(pciv >= 0));                      % 0 or positive values
            assert(all(pciv <= pcOG(value(sgmx))));      % max pc not exceeded
            
        else
            pc = pcOG(s);                      % all cells on drainage curve
        end
    end


end
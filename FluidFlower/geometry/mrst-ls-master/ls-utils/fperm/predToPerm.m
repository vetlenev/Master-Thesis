function k = predToPerm(G, rock, fault, varargin)
% Compute permeability in all directions based on SGR at a particular point
% along a fault.
%
% Input:
% predictor =
% throw     =       
%
% Optional arguments:
%
%
% Output:
% k = permeability across for each fault cell. If k_along is set to true, 
%     then permeability in the y (along horizontal plane) and z (dip) 
%     directions is also given. As a result, k is a double array with size =
%     number of fault cells x 1, or number of fault cells x 3 if k_along 
%     is set to true.
%
%__________________________________________________________________________

% Get variable inputs
opt = struct('name', 'spe02', ...
             'adjust', [], ...
             'kAlong', false, ...
             'kAniso', [], ...
             'cutoff', [], ...
             'smooth', false); 
opt = merge_options(opt,varargin{:});

% Compute permeability
k.name = opt.name;
predictor = fault.kPred;
assert(all(fault.throw > 0));                                               % throw of fault cells must be > 0
switch opt.name
    case 'man99' % Manzocchi et al. (1999)
        % For fault displacement D we use throw, given that the references
        % used by Manzocchi et al. (1999) to get datapoints indicate fault
        % throw whenever they specify what displacement refers to (eg
        % Foxford et al., 1998).
        D = fault.throw;
        SGR = predictor.SGR;
        log_k = -4*SGR - 0.25*log10(D).*(1-SGR).^5;                         % log k [mD]
        k_x = 10.^log_k;                                                    % [mD]
    case 'spe02'
        SGR  = predictor.SGR;
        zmax = G.cells.zmax(fault.fcells);    
        zf   = fault.zf;
        a    = [8*10^4, 19.4, 0.00403, 0.0055, 12.5];
        k_x  = a(1)*exp(-(a(2)*SGR + a(3)*zmax + ...
                          (a(4)*zf - a(5)).*(1-SGR).^7));                   % [mD]
end
k_x = k_x*milli*darcy;

if ~isempty(opt.cutoff)
    switch opt.cutoff
        case 'absMax'
            kAbsMax = max(rock.perm(:, 1));
            k_x(k_x > kAbsMax) = kAbsMax;
        case 'unitMax'
            fc = fault.fcells;
            uul = unique(G.cells.juxt(:, 1)); uur = unique(G.cells.juxt(:, 2));
            uul = uul(uul>0); uur = uur(uur>0);                             %n = 1 is 0 (no fault cell)
            kUnitMax = zeros(G.cells.num, 2);
            for n = 1:numel(uul)                                            % FW                     
                idc = ismember(G.cells.unit(:,1), uul(n));
                idcf = ismember(G.cells.juxt(:,1), uul(n));
                kUnitMax(idcf, 1) = max(rock.perm(idc, 1));
            end
            for n = 1:numel(uur)                                            % HW           
                idc = ismember(G.cells.unit(:,1), uur(n));
                idcf = ismember(G.cells.juxt(:,2), uur(n));
                kUnitMax(idcf, 2) = max(rock.perm(idc, 1));
            end
            kUnitMax(:, 1)  = max(kUnitMax, [], 2);
            kUnitMax(fc, 2) = k_x;
            kUnitMax        = kUnitMax(any(kUnitMax, 2), :);
            k_x             = min(kUnitMax, [], 2);
            % sort k_x from ascending to original order
            [~, id] = sort(fc);
            k_x(id) = k_x;
        case '50'
            kAbsMax = 50*milli*darcy;
            k_x(k_x > kAbsMax) = kAbsMax;
            
        case '200'
            kAbsMax = 200*milli*darcy;
            k_x(k_x > kAbsMax) = kAbsMax;
            
    end
end

if opt.smooth == true
   k_x = smooth(G.cells.centroids(fault.fcells, 3) , k_x, 0.02, 'rloess'); 
end

if opt.kAlong == true
    kmult = opt.kAniso;
    k.vals = [k_x, k_x*kmult, k_x*kmult];                
else
    k.vals = k_x;
end
        
if ~isempty(opt.adjust)
    k.vals = k.vals*opt.adjust;
end

return
end



function [val] = relatePoroPerm(G, rock, fault, varargin)
%
%
%

% Get variable inputs
opt = struct('method', 'carKoz', ...
             'out', 'poro', ...
             'endSand', [], ...
             'endClay', [], ...
             'makePlot', 0); 
opt = merge_options(opt,varargin{:});
latx = {'interpreter', 'latex'};
axfs = {'FontSize', 12};
legfs = {'FontSize', 12};

switch opt.method
    case 'carKoz'
        if strcmp(opt.out, 'poro')
            
        else
            
        end
        error('Carman-Kozeny model not coded yet.')
    
    case 'rev99End'                                 % just for endmembers, test model
        if strcmp(opt.out, 'poro')
            if numel(fault.k.vals(1, :)) > 1
                inProp = geomean(fault.k.vals(:, 1:2), 2);
            else
                inProp = fault.k.vals(:, 1);
            end
            d      = 0.25*10^(-3);                  % sand grain diameter
            m      = [1.7, 3.2];                    % [sand, shale]
            refk   = 0.01*milli*darcy;              % reference shale perm
            refphi = 0.15;                          % reference shale poro
            vcl    = fault.kPred.SGR;               % shale fraction (vol.)
            fc     = fault.fcells;
            
            % Compute
            val = zeros(numel(fc), 1);
            idsand = vcl < 0.7;
            val(idsand) = nthroot(inProp(idsand)*24./(d^2), 3*m(1));
            val(~idsand) = nthroot(inProp(~idsand)./refk, ...
                                   3*m(2)) * refphi;   
            
            % Set cutoffs
            val(val<0.08) = 0.08;
            
            % Add throw dependency (if throw -> 0, poro -> poro juxtaposed
            % unit.
            noThrow = fault.throw < 1;
            cellIds = fc(noThrow);
            uul = unique(G.cells.juxt(cellIds, 1)); 
            for n = 1:numel(uul)                                            % FW 
                idc  = ismember(G.cells.unit(:,1), uul(n));
                idcf = ismember(G.cells.juxt(fc,1), uul(n));
                val(all([noThrow,idcf], 2)) = mean(rock.poro(idc));
            end
            
            % Smooth
            val = smooth(G.cells.centroids(fault.fcells, 3) , val, 0.01, 'rloess');  
        end
        
    case 'rev99'
        
    case 'idealPack'                                % ideal packing. Eg Marion et al. (1992)
        if strcmp(opt.out, 'poro')
            vcl      = fault.kPred.SGR;             % clay fraction (vol.)
            poroSand = opt.endSand(G.cells.centroids(fault.fcells, G.griddim));
            poroClay = opt.endClay(G.cells.centroids(fault.fcells, G.griddim));
            poroSand = reshape(poroSand,numel(fault.fcells), 1);
            poroClay = reshape(poroClay,numel(fault.fcells), 1);
            
            val = zeros(numel(fault.fcells), 1);
            shalySand = vcl < poroSand;
            val(shalySand) = poroSand(shalySand) - vcl(shalySand).*(1-poroClay(shalySand));
            val(~shalySand) = vcl(~shalySand).*poroClay(~shalySand);  
        else
            error('"idealPack" is not a poro-perm relationship. It is to compute porosity only.')
        end
        
        if opt.makePlot == 1
            figure
            scatter(vcl, val, 4, 'k')
            xlabel('$V_\mathrm{cl}$ [-]', axfs{:}, latx{:})
            ylabel(' $\phi$ [-]', axfs{:}, latx{:})
            grid on            
        end
end

end
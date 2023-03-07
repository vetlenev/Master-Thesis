classdef MyPvMult < StateFunction
    % Effective pore-volume after pressure-multiplier
    properties
    end
    
    methods
        function gp = MyPvMult(model, varargin)
            gp@StateFunction(model, varargin{:});
            if isfield(model.fluid, 'pvMultR')
                gp = gp.dependsOn({'pressure'}, 'state');
            end
        end
        function pv = evaluateOnDomain(prop, model, state)
            % Get effective pore-volume, accounting for possible
            % rock-compressibility
            f = model.fluid;
            pv = model.operators.pv;
            p0 = model.operators.p0; %LS inital p (i.e. after initialization)
            if isfield(f, 'pvMultR')
                p  = model.getProp(state, 'pressure');
                %pvMult = prop.evaluateFunctionOnDomainWithArguments(f.pvMultR, p, p0);
                regions = model.rock.regions.rocknum;
                pvMult = evaluateFunctionCellSubsetReg(prop, f.pvMultR, regions, p, p0);
                pv = pv.*pvMult;
            end
        end
    end
end
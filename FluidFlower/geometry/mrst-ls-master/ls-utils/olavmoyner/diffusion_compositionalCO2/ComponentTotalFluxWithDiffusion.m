classdef ComponentTotalFluxWithDiffusion < ComponentTotalFlux
    properties
        D
        componentDiffusion
        faceAverage = false;
    end
    
    methods
        function cf = ComponentTotalFluxWithDiffusion(model)
            cf@ComponentTotalFlux(model);
            ncomp = model.getNumberOfComponents();
            cf.D = getFaceDiffusivity(model.G, model.rock);
            cf.componentDiffusion = ones(1, ncomp);
        end

        function v = evaluateOnDomain(prop, model, state)
            op = model.operators;
            T = prop.D(op.internalConn);
            faceAvg = prop.faceAverage;
            cellMass = model.getProps(state, 'ComponentTotalMass');
            frac = cellMass;
            totMass = 0;
            ncomp = model.getNumberOfComponents();
            for i = 1:ncomp
                totMass = totMass + cellMass{i};
            end
            for i = 1:ncomp
                frac{i} = frac{i}./totMass;
            end

            v = evaluateOnDomain@ComponentTotalFlux(prop, model, state);
            C = prop.componentDiffusion;
            for c = 1:ncomp
                Tc = C(c).*T;
                if ~isempty(v{c})
                    xi = frac{c};
                    grad_xi = op.Grad(xi);
                    if faceAvg
                        faceMass = op.faceAvg(cellMass{c});
                    else
                        flag = value(grad_xi) < 0;
                        faceMass = op.faceUpstr(flag, cellMass{c});
                    end
                    diffuse_flux = -faceMass.*Tc.*grad_xi;
                    v{c} = v{c} + diffuse_flux;
                end
            end
        end
    end
end
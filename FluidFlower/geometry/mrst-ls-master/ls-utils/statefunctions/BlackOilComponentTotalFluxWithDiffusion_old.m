classdef BlackOilComponentTotalFluxWithDiffusion < ComponentTotalFlux
    properties
        D
        componentDiffusion
        faceAverage = false;
    end
    
    methods
        function cf = BlackOilComponentTotalFluxWithDiffusion(model)
            cf@ComponentTotalFlux(model);
            ncomp = model.getNumberOfComponents();
            cf.D = getFaceDiffusivity(model.G, model.rock);
            cf.componentDiffusion = ones(1, ncomp);
        end

        function v = evaluateOnDomain(prop, model, state)
            op = model.operators;
            T = prop.D(op.internalConn);
            faceAvg = prop.faceAverage;
            cellMass = model.getProps(state, 'ComponentPhaseMass');
            frac = cellMass;
            totMass = 0;
            ncomp = model.getNumberOfComponents();
            for i = 1:ncomp
                totMass = totMass + cellMass{i};
            end
            for i = 1:ncomp
                frac{i} = frac{i}./totMass;
            end

            waterInBrine = cellMass{1, 1};
            Co2InBrine = cellMass{2, 1};
            massBrine = waterInBrine + Co2InBrine;

            X_co2 = Co2InBrine./massBrine;
            

            v = evaluateOnDomain@ComponentTotalFlux(prop, model, state);
            C = prop.componentDiffusion;
            oix = model.getPhaseIndex('O');
            cellDens = model.getProps(state, 'Density');
            s = state.s;
            if false
                dix = fliplr(1:ncomp);
                for c = 2%1:ncomp
                    Tc = C(c).*T;
                    if ~isempty(v{c})
                        xi = frac{c};
                        grad_xi = op.Grad(xi);
                        if faceAvg
                            %faceMass = op.faceAvg(cellMass{c});
                            faceDens = op.faceAvg(cellDens{dix(c)});
                            faceS = op.faceAvg(s{dix(c)});
                        else
                            flag = value(grad_xi) < 0;
                            %faceMass = op.faceUpstr(flag, cellMass{c});
                            faceDens = op.faceUpstr(flag, cellDens{dix(c)});
                            faceS = op.faceUpstr(flag, s{dix(c)});
                        end
                        %diffuse_flux = -faceMass.*Tc.*grad_xi;
                        diffuse_flux = -faceDens.*faceS.*Tc.*grad_xi;
                        v{c} = v{c} + diffuse_flux;
                    end
                end
            else
                c = 2;
                grad_xi = op.Grad(X_co2);
                if faceAvg
                    faceDens = op.faceAvg(cellDens{1});
                    faceS = op.faceAvg(s{1});
                else
                    flag = value(grad_xi) < 0;
                    faceDens = op.faceUpstr(flag, cellDens{1});
                    faceS = op.faceUpstr(flag, s{1});
                end
                Tc = C(c).*T;
                diffuse_flux2 = -faceDens.*faceS.*Tc.*grad_xi;
                v{2} = v{2} + diffuse_flux2;
%                 if isa(diffuse_flux2, 'GenericAD')
%                     id = ~isnan(diffuse_flux2.val);
%                     assert(sum(diffuse_flux2.val(id)-diffuse_flux.val(id)) == 0)
%                 else
%                     id = ~isnan(diffuse_flux2);
%                     assert(sum(diffuse_flux2(id)-diffuse_flux(id)) == 0)
%                 end
            end
        end
    end
end
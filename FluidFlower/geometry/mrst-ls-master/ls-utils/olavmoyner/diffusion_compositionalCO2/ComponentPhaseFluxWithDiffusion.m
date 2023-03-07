classdef ComponentPhaseFluxWithDiffusion < ComponentPhaseFlux
    properties
        D
        componentDiffusion
        faceAverage = false;
    end
    
    methods
        function cf = ComponentPhaseFluxWithDiffusion(model)
            cf@ComponentPhaseFlux(model);
            cf = cf.dependsOn({'PermeabilityPotentialGradient', 'FaceComponentMobility'});
            
            ncomp = model.getNumberOfComponents();
            cf.D = getFaceDiffusivity(model.G, model.rock);
            cf.componentDiffusion = ones(1, ncomp);
        end

        function v = evaluateOnDomain(prop, model, state)
            op = model.operators;
            T = prop.D(op.internalConn);
            faceAvg = prop.faceAverage;
            M = model.getProps(state, 'ComponentPhaseMoleFractions');
%             M = model.getProps(state, 'ComponentPhaseDensity');
%             M = model.getProps(state, 'ComponentPhaseMoleFractions');
            cellMass = model.getProps(state, 'ComponentPhaseMass');
            ncomp = model.getNumberOfComponents;
            nph = model.getNumberOfPhases;
            v = evaluateOnDomain@ComponentPhaseFlux(prop, model, state);
            C = prop.componentDiffusion;
            for c = 1:ncomp
                Tc = C(c).*T;
                for ph = 1:nph
                    if ~isempty(v{c, ph})
                        xi = M{c, ph};
                        grad_xi = op.Grad(xi);
                        if faceAvg
                            faceMass = op.faceAvg(cellMass{c, ph});
                        else
                            flag = value(grad_xi) < 0;
                            faceMass = op.faceUpstr(flag, cellMass{c, ph});
                        end
                        diffuse_flux = -faceMass.*Tc.*grad_xi;
                        v{c, ph} = v{c, ph} + diffuse_flux;
                    end
                end
            end
        end
    end
end
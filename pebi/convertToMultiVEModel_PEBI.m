function [model_ve, model_coarse] = convertToMultiVEModel_PEBI(model, varargin)
%Convert a regular model to a hybrid/multilayer VE model
%
% SYNOPSIS:
%   [model_ve, model_coarse] = convertToMultiVEModel(model)
%   [model_ve, model_coarse] = convertToMultiVEModel(model, isFineCells)
%
% REQUIRED PARAMETERS:
%   model  - Fine-scale description in the form of typically a water gas
%            model from the co2lab module.
%   isFine - Optionally an indicator that is true for cells to be retained
%            with a fine discretization in the new model
%
% OPTIONAL PARAMETERS:
%   sealingFaces - List of sealing faces impermeable to flow. These divide
%                  VE regions. If your model has multiple layers you really
%                  should specify this to get good results.
%   
%   sealingCells - List of sealing cells, the (almost) impermeable layer
%                  delimited by sealing faces.
%
%   sealingCells_faces - List of sealing faces delimiting sealing cells.
%                           Are NOT assigned transmissibility multipler,
%                           only used for categorization of hybrid grid.
%
%   setSubColumnsFine  - Boolean indicating if assigning all cells in
%                           subcolumns of imposed fine cells to fine
%                           (true) or only the imposed fine cells (false).
%
%   multiplier   - Weighting applied to sealing faces (default: zero)
%
%   
%   sumTrans     - Upscale coarse transmissibility by summing up values
%                  over each face. Default: true. Otherwse, will use
%                  whatever is the current default in upscaleModel
%
%   transThreshold - Threshold to consider transmissibilities as "sealing"
%                    if sealingFaces is defaulted.
%
%   ve_cols         - cell array of fine cells connected to each VE column
%   
%   pe_rest - Entry pressure in high-permeable region (where VE columns)
%
% RETURNS:
%   model_ve     - A VE model on the coarse scale
% 
%   model_coarse - The same fine-scale discretization, on the coarse scale.
%                  Mostly useful for comparison to the VE version.
% EXAMPLE:
%   introHybridVE
%
% SEE ALSO:
%   convertMultiVEStates

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    if mod(numel(varargin), 2) == 1
        isFine = varargin{1};
        varargin = varargin(2:end);

        if islogical(isFine)
            assert(numel(isFine) == model.G.cells.num, 'Logical mask must match grid size');
            isFine = find(isFine);
        end
    else
        isFine = [];
    end
    opt = struct('remFineFaces', [], 'sealingFaces', [], 'sealingXFaces', [], 'sealingCells', [], 'sealingCells_faces', [], 'setSubColumnsFine', true, ...
                'multiplier', 0, 'sumTrans', true, 'transThreshold', 0, 've_cols', [], 'fine_cols', [], 'pe_rest', 0);
    opt = merge_options(opt, varargin{:});

    G = model.G;
    fluid = model.fluid;
    % Ensure transmissibilities are up-to-date
    trans = model.operators.T_all;
    trans(model.operators.internalConn) = model.operators.T;
    
    % Find categories
    trans_category = ones(model.G.faces.num, 1);
    trans_category(trans <= opt.transThreshold) = 0;
    trans_test = ones(model.G.faces.num, 1);
    trans_test(trans <= opt.transThreshold) = 0;
    % Zero category: different discretization, virtual fine cells
    if ~isempty(opt.sealingFaces) % vertical VE trans
        %trans_category(opt.sealingFaces) = 0;
    end
    % -----
    if ~isempty(opt.sealingCells_faces) % VE to Fine (horizontal or vertical)
        %trans_category(opt.sealingCells_faces) = 0;
    end
    if ~isempty(opt.remFineFaces) % Remaining fine cells (well, bc, ...). NB: transmissibility is NOT set to zero for these cells, only trans_category used to define different discretizations
       trans_category(opt.remFineFaces) = 0;
    end
    % -----

    if ~isempty(opt.sealingXFaces)
        trans_test(opt.sealingXFaces) = 0;
    end

    [categories, c_h, colNo] = findCategoriesMultiVE_PEBI(G, trans_category, opt.ve_cols, opt.fine_cols);
    if ~isempty(isFine) % internal fine -> each fine cell is its own "column"
        if opt.setSubColumnsFine
            % Set category zero for sub-columns containing fine cells
            categories(ismember(c_h, c_h(isFine))) = 0;
        else
            % Only set category zero to imposed fine cells
            categories(isFine) = 0;
        end
    end   
    % Set up grids
    CG = generateCoarseGridMultiVE_PEBI(G, categories, colNo);
    
    %[ii, jj] = gridLogicalIndices(G);
    %colNo =  ii + G.cartDims(1).*(jj-1);
    coarseColNo = zeros(CG.cells.num, 1);
    coarseColNo(CG.partition) = colNo;
    CG.cells.columns = coarseColNo;

    % Extrude to 3D
    if CG.griddim == 2
        CG.cells.centroids(:,3) = 0;
        CG.faces.centroids(:,3) = 0;
        %CG.nodes.coords(:,3) = 0;
    end
    
    disp('before sumTrans')
    disp(min(trans))
      
    if opt.sumTrans
        % Upscale -> average transmissibilities        
        faceno = rldecode(1:CG.faces.num, diff(CG.faces.connPos), 2)';   
        T_c = accumarray(faceno, trans(CG.faces.fconn));      
        %T_c = 1./accumarray(faceno, 1./trans(CG.faces.fconn));
        
        % Generate a coarse model 
        model_coarse = upscaleModelTPFA_PEBI(model, CG.partition, 'transCoarse', T_c);
    else     
        model_coarse = upscaleModelTPFA_PEBI(model, CG.partition); % use coarse rock to upscale trans
    end
    
    disp('after sumTrans')
    disp(min(model_coarse.operators.T_all))
    % Grab rock from coarse model
    rock_c = model_coarse.rock;
     
    % Copy operators over and create a new VE model
    if isa(model, 'TwoPhaseWaterGasModel_PEBI')
        model_ve = WaterGasMultiVEModel_PEBI(CG, rock_c, fluid, 'sealingFaces', opt.sealingFaces, 'sealingCells', opt.sealingCells, 'pe_rest', opt.pe_rest);
    elseif isa(model, 'OverallCompositionCompositionalModel')
        model_ve = OverallCompositionMultiVEModel(CG, rock_c, fluid, model.EOSModel);
    elseif isa(model, 'NaturalVariablesCompositionalModel')
        model_ve = NaturalVariablesMultiVEModel(CG, rock_c, fluid, model.EOSModel);
    else
        error(['VE not implemented for class ', class(model)]);
    end
    badc = model_ve.G.cells.height <= 0;
    mv = max(1e-4*mean(model_ve.G.cells.height), 1e-10);
    if any(badc)
        warning(['Found ', num2str(nnz(badc)), ' cells with zero height. Fixing...']);
        model_ve.G.cells.height(badc) = mv;
        model_ve.G.cells.bottomDepth(badc) = model_ve.G.cells.bottomDepth(badc) + mv;
    end
    delta = model_ve.operators.connections.faceBottomDepth - model_ve.operators.connections.faceTopDepth;
    badf = model_ve.operators.connections.faceHeight(:) <= 0 | delta(:) <= 0;
    if any(badf)
        warning(['Found ', num2str(nnz(badf)), ' faces with zero height. Fixing...']);
        model_ve.operators.connections.faceHeight(badf) = mv;
        model_ve.operators.connections.faceBottomDepth(badf) = model_ve.operators.connections.faceBottomDepth(badf) + mv;
    end
       
    
    if ~isempty(opt.sealingFaces) % && ~opt.sumTrans
        coarseFaceNo = rldecode((1:CG.faces.num)', diff(CG.faces.connPos));
        toMult = false(CG.faces.num, 1);
        if islogical(opt.sealingFaces)
            map = opt.sealingFaces;
        else
            map = false(G.faces.num, 1);
            map(opt.sealingFaces) = true;
        end
        
        toMult(coarseFaceNo(map(CG.faces.fconn))) = true;
        model_coarse.operators.T_all(toMult) = model_coarse.operators.T_all(toMult).*opt.multiplier;
        model_coarse.operators.T = model_coarse.operators.T_all(model_coarse.operators.internalConn);
    end     

    model_ve.nonlinearTolerance = model.nonlinearTolerance;
    model_ve.minimumPressure = model.minimumPressure;
    model_ve.extraStateOutput = model.extraStateOutput;

    model_ve.operators.T = model_coarse.operators.T;
    model_ve.operators.T_all = model_coarse.operators.T_all;
    
    % --- Add veBottom cells and connections ---
    op = model_ve.operators;
    isVE = model_ve.G.cells.discretization > 1;
    veTransition = op.connections.veToFineVertical | ...
                    op.connections.veTransitionVerticalConn & op.T > 0;        
    veAll = op.connections.veInternalConn | op.connections.veTransitionHorizontalConn;
    cn = op.N(veTransition, :);
    c_vic = op.N(veAll, :);  
   
    cB = [];
    veB = [];
    for idx=1:2
        c = cn(:,idx);
        isVE_c = isVE(c);
        if any(isVE_c)
            t = model_ve.G.cells.topDepth(c);
            T = op.connections.faceTopDepth(veTransition, idx);
            b = model_ve.G.cells.bottomDepth(c);
            B = op.connections.faceBottomDepth(veTransition, idx);
            Hb = model_ve.G.cells.height(c); % to not overwrite global H
            
            veB = B == b & T ~= t;
            cB = c(veB); % select correct bottom cells
%             cb = c_vic(ismember(c_vic, cb));
%             if ~isempty(cb)
%                 cB = cat(2, cB, cb);
%             end
        end
    end    
    
    model_ve.G.cells.bottomVE = cB;
    model_ve.operators.connections.bottomVE = veB;
    % ------------------------------------------
end
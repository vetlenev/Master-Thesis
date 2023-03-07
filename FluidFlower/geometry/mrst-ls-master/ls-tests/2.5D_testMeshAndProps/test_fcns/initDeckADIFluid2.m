function fluid = initDeckADIFluid2(deck, varargin)
%Initialize AD-solver fluid from ECLIPSE-style input deck
%
% SYNOPSIS:
%   f = initDeckADIFluid(deck)
%
% REQUIRED PARAMETERS:
%   deck - Output structure from readEclipseDeck
%
% OPTIONAL PARAMETERS:
%   G    - Grid to be used for the fluid model. Only required for the
%          situation when a deck has more than one region, and the intended
%          simulation grid is different from the one defined by the deck
%          (more specifically, the deck.GRID.ACTNUM field).
%
%   region_method - Method for defining regions. Either 'deck' (use exactly as
%           prescribed in the input deck [DEFAULT] or 'override' which
%           allows individual fields to be overwritten. Since the region
%           support in MRST is abstracted away, this option should only be
%           used if you are comfortable with the internal structure of
%           getRegMap and interpReg.
%
%   regionOverride - Struct containing fields to be overwritten. Only used
%           if 'method' equals 'override'. The supported fields include
%           PVTNUM, SATNUM, ROCKNUM, IMBNUM and SURFNUM. For each of these
%           fields, if present, the value can be formatted in three
%           different ways:
%               * Single numerical value. Will be repeated for all grid
%                 cells.
%               * The char ':'. Interpreted as using the first region for
%                 all cells. This is normally used when no regions are
%                 present.
%               * A Nx1 array containing one region indicator per cell in
%                 the grid. The region indicator is a single number
%                 indicating which table is to be used for that cell.
%
%   singleRegion - Internal debug option. This region will be used for ALL
%                 fields and any other value than 1 or ':' will usually
%                 result in exceptions being thrown.
%
% RETURNS:
%   f - Fluid model suitable for several different MRST models.
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, TwoPhaseOilWaterModel, initEclipseDeck

%{
Copyright 2009-2018 SINTEF Digital, Mathematics & Cybernetics.

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

% this only work for full deck or first region.
opt = struct('G',              [], ...
             'region_method', 'deck', ...
             'singleRegion',   [], ...
             'regionOverride', struct());
 
opt = merge_options(opt,varargin{:});
% reg = handleRegions(deck, opt.G);
reg = getRegions(deck, opt);

fluid = struct();

% Properties
props = deck.PROPS;

for fld = assignable_fields(fieldnames(props))
   try
      fluid = feval(['assign', fld{1}], fluid, props.(fld{1}), reg);

   catch ME
      warning(msgid('Assign:Failed'), ...
             ['Could not assign property ''%s''. ', ...
              'Encountered error: ''%s'''], fld{1}, ME.message);
   end
end

fn = fieldnames(fluid);
for i = 1:numel(fn)
    f = fn{i};
    if iscell(fluid.(f)) && numel(fluid.(f)) == 1
        fluid.(f) = fluid.(f){1};
    end
end
end

function reg = getRegions(deck, opt)
    % Legacy regions
    reg = handleRegions(deck, opt.G);
    % Modern region treatment
    reg.sat = 1;
    reg.pvt = 1;
    reg.rock = 1;
    if isfield(deck.RUNSPEC, 'TABDIMS')
        tab = deck.RUNSPEC.TABDIMS;
        reg.sat = tab(1);
        reg.pvt = tab(2);
        reg.rock = tab(7);
    end
end

%--------------------------------------------------------------------------

function flds = assignable_fields(fnames)
   excpt = reshape(exclude_properties(), [], 1);
   n     = numel(excpt);

   [i, j] = blockDiagIndex(numel(fnames), n);

   mtch = strcmp(reshape(fnames(i), [], n), ...
                 reshape(excpt (j), [], n));

   flds = reshape(fnames(~ any(mtch, 2)), 1, []);
end

%--------------------------------------------------------------------------

function excpt = exclude_properties()
   % Properties not resulting in individual functions
   excpt = {'KRW'   , 'KRO'   , 'KRG',   ...
            'SWL'   , 'SWCR'  , 'SWU',   ...
            'SGL'   , 'SGCR'  , 'SGU',   ...
            'SOWCR' , 'SOGCR' ,          ...
            'ISWL'  , 'ISWCR' , 'ISWU',  ...
            'ISGL'  , 'ISGCR' , 'ISGU',  ...
            'ISOWCR', 'ISOGCR',          ...
            'CNAMES', 'BIC'   , 'ACF',   ...
            'PCRIT' , 'TCRIT' , 'VCRIT', ...
            'MW'    , 'ZCRIT' , 'EOS'    ...
            'SCALECRS', 'SWATINIT', ...
            };
end

%--------------------------------------------------------------------------

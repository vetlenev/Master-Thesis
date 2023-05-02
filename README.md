# Master-Thesis
Student: Vetle Nevland.

Supervisors: Knut-Andreas Lie, Odd Andersen.

Codes and files relevant for my master thesis on hybrid modelling of CO2 storage with applications to FluidFlower (https://fluidflower.w.uib.no/).

Folders:

`FluidFlower`
- Grid generation and hybrid modelling of FluidFlower.
- Geometric data of layering of FluidFlower rig provided in `geometry` folder.
- A semi-structured grid for FluidFlower is generated by the classes and functions in `grid` folder.
- `discretization` and `models` includes governing equations, discretizations and utilities for setting up a hybrid model compatible with the semi-structured grid.

`hybrid2D`
- Codes for setting up simulations on synthetic, two-dimensional heterogeneous formations.

`hybrid3D`
- Codes for setting up simulations on synthetic, three-dimensional heterogeneous formations.

`stochastic`
- Folder to store seeds and figures generated for simulations on two-dimensional heterogeneous formation with stochastic distribution of semi-permeable layers. 

`trapping`
- Codes extending the trap analysis from [mrst-co2lab](https://www.sintef.no/projectweb/mrst/modules/co2lab/) to be applicable for hybrid modelling.
- Also includes trap analysis for full-dimensional model for comparison of models.

`models`
- Extended hybrid model implementing residual saturation and capillary exclusion in a relaxed VE setting.

`discretization`
- Discretization of governing equations on coarse and fine scales, including treatment of transitions between different discretization regions.

Link to thesis project:
- Under construction ...

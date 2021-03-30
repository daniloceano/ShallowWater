# ShallowWater

2D Limited-area shallow-water model with nonlinear advection terms.

The program was developed as an exercise during the course on Atmospherical Modelling, at the Department of Astronomy, Geosciences and Atmospherical Sciences of The University of São Paulo (USP), ministrated by professor Pedro Leite da Silva Dias, 2020.

The references for the formulas can be found at: Basic Numerical Methods in Meteorology and Oceanography, by Kristofer Döös, Peter Lundberg and Aitor Aldama Campino.


The program core.py calls the subtourines:

* grid.py, for creating the computational grid
* forcing.py, for setting the initial conditions, in this case a forcing in the meridional wind component that slowly peaks and decays
* diagnostic_eqs.py, for calculating the mass, velocty and vorticity fluxes
* filters.py, for applying filters in order to dump computational modes that appear during model intergation
* 4 distinct boundaries routines, that use distinct methods for calculating radation boundary conditions, i. e., allowing the solution to leave the computational domain
*  ex5_plots.py, which  stores the model solution for a set of experiments and makes figures
*  filme.py, which takes the outputs from ex5_plots.py and create animations

# TASKit
## A Python-based kit containing Post-MD analysis tools.

Current features are: 
* Generation of potential energy timeseries from MD, using OpenMM
* Trajectory decorrelation, using J. Chodera's PyMBAR library
* RMSD All-vs-All computation and generation of corresponding heatmap, using MDTraj
* Clusterization (and visualisation using Plotly) and extraction of representative structures, using SciPy 


## Installation

Since it's just a series of python modules, there is no compilation needed. However, the location of the "TASKit" folder needs to be found by your python interpreter. To do that, either add it to your $PYTHONPATH environment variable or to your sys.path variable inside of python.

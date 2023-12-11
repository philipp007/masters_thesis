Masters thesis project 'Modelling and coupled simulation of a thermal and electrochemical system' by Philipp Kretz as part of the program 'Computational Science and Engineering' at TUM. 

This repository contains all software written in cooperation with the company Knorr-Bremse AG and as part of one of their projects. 
Please contact me in case of questions: philipp.kretz@tum.de


Required software: 
- Ubuntu 22.04 or newer
- Python 3 or newer
- preCICE 2.5.0 or newer ([link](https://precice.org/installation-packages.html))
- legacy FEniCS 2019.1.0 ([link](https://fenicsproject.org/download/archive/))
- meshio ([link](https://pypi.org/project/meshio/))
- (Gmsh 4.11.1 ([link](https://gmsh.info/doc/texinfo/gmsh.html)))
- (paraview ([link](https://www.paraview.org/download/)))


Built: 
there is no built required 

Run the software: 
measurement_data/pp_dummy.py:
- plug in path to .csv-file which should be processed in line 304
- add any desired output to end of the script (figures showing the original data and processed data are shown automatically)
- starte the script by calling 'python3 ./pp_dummy.py'

battery simulation: 
- modify battery properties in thermal/thermal_model.py, electrostatics/electrostatic_model.py and utilities/battery_properties.py (class Battery)
- modify physical domains for simulations using Gmsh
  - electrostatics:  import 25d_mesh.geo and export 25d_mesh.msh-file
  - thermal:         import toshiba_cell.geo and export toshiba_cell.msh-file
- open two terminal windows and cd to folders thermal and electrostatics respectively
  - call 'python3 ./thermal_model.py' in the first window
  - call 'python3 ./electrostatic_model.py' in the second window

Ouput of the battery simulation: 
- graph of battery properties pops up automatically
- output-files of electrostatic and thermal simulation are stored in the respective 'output'-folder and can be visualized e.g. by using paraview
- there is no post-processing needed

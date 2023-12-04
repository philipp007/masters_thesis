from fenics import *
import meshio 
import numpy as np
import sys
sys.path.append('../utilities')
from battery_properties import *
import precice
#from fenicsprecice import Adapter
import matplotlib.pyplot as plt


## parameters
fenics_dt = 10                  #s
dt = fenics_dt
T_amb = 293.15                  #K
Heat_cap = 500                  #J/K

tc = [25, 25, 2]             #W/m/K


class Thermal_Conductivity(UserExpression):
    def __init__(self, tc, **kwargs):
        super().__init__(**kwargs)
        self.tc = tc

    def eval_cell(self, values, x, cell):
        k = np.array([[self.tc[0], 0, 0],
                      [0, self.tc[1], 0],
                      [0, 0, self.tc[2]]
                      ])
        values[:] = k.flatten()

    def value_shape(self):
        return(3, 3)



thermal_conductivity = Thermal_Conductivity(tc)


## mesh import 
msh = meshio.read("toshiba_cell.msh")
meshio.write("dolfin_files/mesh.xdmf", meshio.Mesh(points=msh.points, \
                                    cells={"tetra": msh.cells["tetra"]}))
Surfaces = ["terminal_negative", 
            "terminal_positive", 
            "housing_I", 
            "housing_VI", 
            "housing_II", 
            "housing_V", 
            "housing_III", 
            "housing_IV"]
for i in range(len(Surfaces)):
    meshio.write("dolfin_files/" + Surfaces[i] + ".xdmf", \
                meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}, \
                cell_data={"triangle": {Surfaces[i]: msh.cell_data["triangle"]["gmsh:physical"]}}))
    
mesh = Mesh()
with XDMFFile("dolfin_files/mesh.xdmf") as infile:
    infile.read(mesh)
mvc = []
meshio.write("dolfin_files/housing_I.xdmf", \
            meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]},
            cell_data={"triangle": {"housing_I": msh.cell_data["triangle"]["gmsh:physical"]}}))

groups = []
for i in range(len(Surfaces)):
    mvc.append(MeshValueCollection("size_t", mesh, 2))
    with XDMFFile("dolfin_files/" + Surfaces[i] + ".xdmf") as infile:
        infile.read(mvc[-1], Surfaces[i])
    tmp_group = cpp.mesh.MeshFunctionSizet(mesh, mvc[-1])
    groups.append(tmp_group)
File("dolfin_files/dolfinmesh.pvd").write(mesh)
boundary_conditions = {}
for i in range(len(Surfaces)):
    File("dolfin_files/dolfincellfunc.pvd").write(groups[i])


mv = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("dolfin_files/housing_I.xdmf") as infile: 
    infile.read(mv, "housing_I")
mf = cpp.mesh.MeshFunctionSizet(mesh, mv)
File("dolfin_files/dolfincellfunc.pvd").write(mf)


## initialize electro thermal interface
eti = electro_thermal_interface()
eti.from_mesh(mesh)

## define Functionspace
V = FunctionSpace(mesh, 'P', 1)

## define Vectorspace for current vector field 
V_vec = VectorFunctionSpace(mesh, "CG", 1)

## define Trial and Test Functions
u = TrialFunction(V)
v = TestFunction(V)
Source = Function(V)
rho_c = Heat_cap / eti.geom_battery.V
u_n = Function(V)

## define boundary conditions
boundary_conditions = \
    {0: {'Dirichlet': T_amb + 80, 'num': 1},      # neg terminal
     1: {'Dirichlet': T_amb + 80, 'num': 2},      # pos terminal
     2: {'Dirichlet': T_amb, 'num': 3},           # housing_I     front 
     3: {'Dirichlet': T_amb, 'num': 4},           # housing_VI    back  
     4: {'Dirichlet': T_amb, 'num': 5},           # housing_II    left
     5: {'Dirichlet': T_amb, 'num': 6},           # housing_V     right
     6: {'Dirichlet': T_amb, 'num': 7},           # housing_III   top 
     7: {'Dirichlet': T_amb, 'num': 8},           # housing_IV    bottom 
     }

bcs = []
ds_C = []
integrals_N = []
for i in boundary_conditions:
    if 'Dirichlet' in boundary_conditions[i]:
        bc = DirichletBC(V, boundary_conditions[i]['Dirichlet'], \
                         groups[i], boundary_conditions[i]['num'])
        
        bcs.append(bc)
    elif 'Neumann' in boundary_conditions[i]:
        if boundary_conditions[i]['Neumann'] != 0:
            g = boundary_conditions[i]['Neumann']
            ds_C.append(Measure("ds", domain=mesh, \
                subdomain_data=groups[i], subdomain_id=boundary_conditions[i]['num']))
            integrals_N.append(g*v*ds(i))
    else:
        raise Exception("unknwon boundary type")

ds_C = (Measure("ds", domain=mesh, subdomain_data=mf, subdomain_id=3))


bcs.append(bc)


## preCICE setup
participant_name = "thermal"
config_file_name = "../precice-config.xml"
solver_process_index = 0
solver_process_size = 1
interface = precice.Interface(participant_name, config_file_name, \
                              solver_process_index, solver_process_size)

mesh_name = "thermal-Mesh"
mesh_id = interface.get_mesh_id(mesh_name)

data_name = "heat"
data_id = interface.get_data_id(data_name, mesh_id)

positions = eti.get_positions()

vertex_ids = interface.set_mesh_vertices(mesh_id, positions)

precice_dt = interface.initialize()


## define Input and weak form 
F =  u_n*v*dx - dt*dot(thermal_conductivity*grad(u), grad(v))*dx \
    - u*v*dx + Source*v*dt/rho_c*dx + sum(integrals_N)
a, L = lhs(F), rhs(F)

## initialize variables
u = Function(V)
u_init = Constant(T_amb)
u_n.assign(u_init)
t = 0
dt = fenics_dt

## initialize output files
vtkfile = File('output/output.pvd')



## main loop 
while interface.is_coupling_ongoing():

    dt = np.min([fenics_dt, precice_dt])
    eti.read_data(interface.read_block_scalar_data(data_id, vertex_ids))
    Source.assign(eti.map_to_fenics(mesh, V))

    ## solve variational problem 
    t += dt
    solve(a == L, u, bcs)
    u_n.assign(u)
    #print(u.vector()[:])
    precice_dt = interface.advance(dt)

    ## output of fields to vtk 
    vtkfile << (u, t)


# Hold plot
interface.finalize()




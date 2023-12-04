from fenics import *
import meshio 
import matplotlib.pyplot as plt
import sys
sys.path.append('../utilities')
from utilities_electrostatics import *
from battery_properties import *
import precice
import numpy as np

def surf_integral_out(I):
    res = 15
    points_1D_y = np.linspace(33, 47, num=res)
    points_1D_z = np.array([0, 0.5])
    
    values = np.zeros(shape=(res, 2))
    for i in range(res):
        for j in range(2):
            [values[i, j], tmp1, tmp2] = I(51.465, points_1D_y[i], points_1D_z[j])

    integral = 0
    dy = 1
    dz = 0.5
    for i in range(res - 1):
        grid_values = [values[i, 0], values[i, 1], \
                       values[i + 1, 0], values[i + 1, 1]]
        integral += interpolate_2D(dy/2, dz/2, dy, dz, grid_values)

    dA = dy*dz
    
    return integral*dA
    
def surf_integral_in(I):
    res = 100
    points_1D_x = np.linspace(-102.73/2, 102.73/2, num=res)
    points_1D_y = np.linspace(-115/2, 115/2, num=res)
    
    values = np.zeros(shape=(res, res))
    for i in range(res):
        for j in range(res):
            [tmp1, tmp2, values[i, j]] = I(points_1D_x[i], points_1D_y[j], 1)

    integral = 0
    dx = 102.73/res
    dy = 115/res
    for i in range(res - 1):
        for j in range(res - 1):
            grid_values = [values[i, j], values[i, j + 1], \
                           values[i + 1, j], values[i + 1, j + 1]]
            integral += interpolate_2D(dx/2, dy/2, dx, dy, grid_values)

    dA = dx*dy
    
    return integral*dA




## varibles 

fenics_dt = 5                       #s
u_in = 0                            #V
I_in = 150                          #A
rho = 20*10**-3                   #ohm*mm
eti = electro_thermal_interface()
battery = Battery(1, 293)           # Wh, K
dbc = MyExpression0(battery)


## mesh import 
msh = meshio.read("25d_mesh.msh")
meshio.write("dolfin_files/mesh.xdmf", meshio.Mesh(points=msh.points, cells={"tetra": msh.cells["tetra"]}))
meshio.write("dolfin_files/mf.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}, 
                                    cell_data={"triangle": {"electrolyte": msh.cell_data["triangle"]["gmsh:physical"]}}))
meshio.write("dolfin_files/mf2.xdmf", meshio.Mesh(points=msh.points, cells={"triangle": msh.cells["triangle"]}, 
                                    cell_data={"triangle": {"terminal": msh.cell_data["triangle"]["gmsh:physical"]}}))
mesh = Mesh()
with XDMFFile("dolfin_files/mesh.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("dolfin_files/mf.xdmf") as infile: 
    infile.read(mvc, "electrolyte")
mf = cpp.mesh.MeshFunctionSizet(mesh, mvc)
with XDMFFile("dolfin_files/mf2.xdmf") as infile: 
    infile.read(mvc, "terminal")
mf2 = cpp.mesh.MeshFunctionSizet(mesh, mvc)
File("dolfin_files/dolfinmesh.pvd").write(mesh)
File("dolfin_files/dolfincellfunc.pvd").write(mf)
File("dolfin_files/dolfincellfunc.pvd").write(mf2)


## initialization of model interface and battery properties
heat = scalar_field()
heat.from_eti(eti)
I_vf = vector_field(eti)
I_sf = scalar_field()
I_sf.from_eti(eti)


## define Functionspace
V = FunctionSpace(mesh, 'P', 1)

## define Vectorspace for current vector field 
V_vec = VectorFunctionSpace(mesh, "CG", 1)


## define boundary conditions
boundary_conditions = {0: {'Dirichlet': dbc, 'group': 45},      #terminal
                       1: {'Neumann': I_in, 'group': 46}}       #interaction surface
bcs = []
bc = DirichletBC(V, boundary_conditions[0]['Dirichlet'], mf2, boundary_conditions[0]['group'])
bcs.append(bc)
ds_C = (Measure("ds", domain=mesh, subdomain_data=mf, subdomain_id=boundary_conditions[1]['group']))


## define Trial and Test Functions
u = TrialFunction(V)
v = TestFunction(V)
u_n = Function(V)


# preCICE setup
participant_name = "electrostatics"
config_file_name = "../precice-config.xml"
solver_process_index = 0
solver_process_size = 1
interface = precice.Interface(participant_name, config_file_name, solver_process_index, solver_process_size)

mesh_name = "electrostatics-Mesh"
mesh_id = interface.get_mesh_id(mesh_name)

data_name = "heat"
data_id = interface.get_data_id(data_name, mesh_id)

positions = eti.get_positions()

vertex_ids  = interface.set_mesh_vertices(mesh_id, positions)

precice_dt = interface.initialize()


## define Input and weak form 
i_in = Constant(I_in/eti.geom_battery.A)
dt = fenics_dt
F = i_in*v*ds - dot(grad(u), grad(v))/rho*dx - u*v/dt*dx + u_n*v/dt*dx
a, L = lhs(F), rhs(F)


## initialize variables
u = Function(V)
u_init = Constant(battery.get_voltage())
u_n.assign(u_init)
t = 0

## initialize output files
vtkfile = File('output/output.pvd')
vtkfile_I = File('output/output_I.pvd')

## prepare output arrays 
out_voltage = [battery.get_voltage()]
out_current = [I_in]
out_SOC = [battery.SOC]
out_heat = [0]
out_time = [0]

## main loop 
while interface.is_coupling_ongoing():

    if battery.SOC > 0.01:      # in case battery is not empty
        ## solve variational problem 
        dt = np.min([fenics_dt, precice_dt])
        t += dt
        solve(a == L, u, bcs)
        u_n.assign(u)


        ## prepare current density for coupling 
        I = project(grad(u), V_vec)
        I_vf.from_fenics(I)
        I_sf.from_vector_field(I_vf)
        Idens = I_sf.scalar_cell_current()

        ## calculate heat field proportional to current density
        heat = Idens * (battery.get_heat(I_in) / Idens.sum())

        ## update battery
        battery.update_SOC(I_in, dt)

    else:                       # in case battery is empty
        heat.data[:, :] = 0
        print("Battery dead")
        I_in = 0

    ## communicate heat to thermal model
    interface.write_block_scalar_data(data_id, vertex_ids, heat.data.flatten())

    ## tell preCICE to update time
    precice_dt = interface.advance(dt)

    ## log battery parameters 
    out_voltage.append(battery.get_voltage())
    out_current.append(I_in)
    out_SOC.append(battery.SOC)
    out_heat.append(battery.get_heat(I_in))
    out_time.append(out_time[-1] + dt)


    ## output of fields to vtk 
    vtkfile << (u, t)
    vtkfile_I << (I, t)


out = output_battery(out_voltage, out_current, out_SOC, out_heat, out_time)
out.plot_me()

# Hold plot
interface.finalize()
import numpy as np
from fenics import *
import sys
import warnings
import matplotlib.pyplot as plt
sys.path.append('../measurement_data/data')
#from utilities_electrostatics import *

class electro_thermal_interface:

    def __init__(self):
        self.num_x = 20
        self.num_y = 20
        self.geom_battery = geom_Battery()
        self.set_positions()


    def from_mesh(self, mesh):
        
        self.coords = mesh.coordinates()
        self.num_vertices = np.zeros(shape=(self.num_x-1, self.num_y-1))
        for v in vertices(mesh):
            v_idx = v.index()
            [i, j] = self.find_cell(self.coords[v_idx])
            self.num_vertices[i, j] += 1
        


    def read_data(self, data):
        self.data = data

    
    def h(self):
        return self.geom_battery.height
    
    def w(self):
        return self.geom_battery.width
    
    def d(self):
        return self.geom_battery.depth
    
    def set_positions(self):
        x = np.linspace(-self.h()/2, self.h()/2, self.num_x)
        y = np.linspace(-self.w()/2, self.w()/2, self.num_y)

        self.positions = [[x0, y0] for x0 in x for y0 in y]

    
    def get_positions(self):
        return self.positions
    
    def get_mesh_info(self):
        return [self.num_x, self.num_y, self.geom_battery.height, self.geom_battery.width]
    

    def map_to_fenics(self, mesh, V):
        
        Source = Function(V)
        
        
        shifted_values = np.zeros_like(self.num_vertices)

        cell_mid_x = self.h() / self.num_x / 2
        cell_mid_y = self.w() / self.num_y / 2

        grid_shape = self.num_vertices.shape
        for i in range(grid_shape[0]):
            for j in range(grid_shape[1]):
                cell_idx = i*self.num_y + j
                tmp_value = self.interpolate(cell_idx, cell_mid_x, cell_mid_y)
                shifted_values[i, j] = tmp_value
                
        sum_shifted_values = sum(sum(shifted_values))

        if sum_shifted_values > 0:
            sum_self_data = sum(self.data)

            for v in vertices(mesh):
                v_idx = v.index()
                [i, j] = self.find_cell(self.coords[v_idx])
                tmp_value = shifted_values[i, j] / sum_shifted_values * sum_self_data / self.num_vertices[i, j]
                Source.vector()[v_idx] = tmp_value

        else: 
            Source = Constant(0)
        
        return Source
    

    def find_cell(self, coords):

        i, j = 1, 1

        while i < self.num_x-1 and coords[0] > self.positions[i*self.num_y][0]:
            i += 1
        i -= 1
        while j < self.num_y-1 and coords[1] > self.positions[i*self.num_y + j][1]:
            j += 1
        j -= 1
        #cell_idx = i*self.num_y
        

        #cell_x = coords[0] - self.positions[cell_idx][0]
        #cell_y = coords[1] - self.positions[cell_idx][1]
        
        return [i, j]


    def interpolate(self, cell_idx, cell_x, cell_y):

        dx = self.h() / self.num_x
        dy = self.w() / self.num_y

        if abs(cell_x) > dx or abs(cell_y) > dy:
            value = 0
        else:
            grid_values = [self.data[cell_idx], \
                self.data[cell_idx + 1], \
                self.data[cell_idx + self.num_y], \
                self.data[cell_idx + self.num_y + 1]]

            value = interpolate_2D(cell_x, cell_y, dx, dy, grid_values)
            
        return value 


def interpolate_2D(cell_x, cell_y, dx, dy, grid_values):
    [val_1, val_2, val_3, val_4] = grid_values

    val_12 = (val_2 - val_1) / dy * cell_y + val_1
    val_34 = (val_4 - val_3) / dy * cell_y + val_3

    value = (val_34 - val_12) / dx * cell_x + val_12 
    
    return value


class geom_Battery:
    def __init__(self):
        self.height = 102.73
        self.width = 115
        self.depth = 21
        self.A = self.height * self.width
        self.V = self.A * self.depth
        self.housing_surf = 2 * self.height * self.width + \
                            2 * self.height * self.depth + \
                            1 * self.width * self.depth



class Battery: 
    def __init__(self, SOC, T_amb):
        self.SOC = SOC
        self.nom_W = 40         #Wh
        self.nom_Cap = 15       #Ah
        self.T_amb = T_amb - 273.15 #conversion from K to ˚C    
        

        self.x_temperature = np.load('../measurement_data/data/x_temperature.npy', allow_pickle=True)
        self.y_C_rate = np.load('../measurement_data/data/y_C_rate.npy', allow_pickle=True)
        self.z_charge_caps = np.load('../measurement_data/data/z_charge_caps.npy', allow_pickle=True)
        self.z_relative_losses = np.load('../measurement_data/data/z_relative_losses.npy', allow_pickle=True)
        self.z_voltage = np.load('../measurement_data/data/z_voltage.npy', allow_pickle=True)
        #self.z_SOC = np.load('../measurement_data/data/z_SOC.npy', allow_pickle=True)
        
        self.calc_voltage(0)


    def grid_cell_info(self, x, y):

        i, j = -1, -1
        data_points = [2, 3, 4, 4]

        for idx in range(len(self.x_temperature)):
            if self.x_temperature[idx] < x:
                i = idx

        for jdx in range(len(self.y_C_rate)):
            if self.y_C_rate[jdx] < y:
                j = jdx

        warn_str = "the model is neither validated nor is extrapolation reasonable, because the battery cell couldn't perform in this scenario"
        if i == -1 or i == 3:
            if x < 50:
                i = 2
                warnings.warn("cell values are being extrapolated reasonably due to high temperature")
            else:
                raise Exception("cell temperature out of bounds, " + warn_str)
        if j == -1 or j >= data_points[i]: 
            if y >= 0:
                j = 0
                warnings.warn("cell values are being extrapolated reasonably due to low C-rate")
            elif y <= 1.2 * self.y_C_rate[data_points[i]]:
                j = data_points[i] - 1
            else:
                raise Exception("cell current/C-rate is out of bounds, " + warn_str)
        

        cell_x = x - self.x_temperature[i]
        cell_y = y - self.y_C_rate[j]
        dx = self.x_temperature[i+1] - self.x_temperature[i]
        dy = self.y_C_rate[j+1] - self.y_C_rate[j]

        return [i, j, cell_x, cell_y, dx, dy]



    def get_relative_losses(self, I_in):
        C_rate = I_in / self.nom_Cap
        [i ,j, cell_x, cell_y, dx, dy] = self.grid_cell_info(self.T_amb, C_rate)

        grid_values = [self.z_relative_losses[i, j], \
                       self.z_relative_losses[i, j+1], \
                       self.z_relative_losses[i+1, j], \
                       self.z_relative_losses[i+1, j+1]]
        
        return interpolate_2D(cell_x, cell_y, dx, dy, grid_values)

    def calc_voltage(self, C_rate):
        [i ,j, cell_x, cell_y, dx, dy] = self.grid_cell_info(self.T_amb, C_rate)
        idx = int(round((1 - self.SOC) * 100))
        grid_values = [self.z_voltage[i, j, idx], \
                       self.z_voltage[i, j+1, idx], \
                       self.z_voltage[i+1, j, idx], \
                       self.z_voltage[i+1, j+1, idx]]

        self.voltage = interpolate_2D(cell_x, cell_y, dx, dy, grid_values)


    def update_SOC(self, I, dt):
        C_rate = I/self.nom_Cap
        self.calc_voltage(C_rate)
        W = self.voltage * I * dt * 1e-6/3.6 #mV*A*s*1e-6/3.6 = Wh
        self.SOC -= W/self.nom_W
        if self.SOC < 0.01:
            warnings.warn("Battery dead: SOC <= 0")


    def get_voltage(self):
        return self.voltage*1e-3
    

    def get_heat(self, I_in):
        # current to heat transfer coefficient 
        P = self.voltage * I_in * 1e-3   #W
        heat = P * self.get_relative_losses(I_in)

        return heat
    

class MyExpression0(UserExpression):
    def __init__(self, battery, **kwargs):
        super().__init__(kwargs)
        self.battery = battery

    def eval(self, value, ufc_cell):
        value[0] = self.battery.get_voltage() / 2

    def value_shape(self):
        return ()
    

class output_battery:

    def __init__ (self, voltage, current, SOC, heat, time):
        self.voltage = voltage
        self.current = current 
        self.SOC = SOC
        self.heat = heat
        self.time = time 



    def my_subplot(self, ax, value, title):
        ax.set_title(title)
        ax.plot(self.time, value)
        ax.locator_params(nbins=3)
        ax.set_xlabel('time / s')
        #max_val = abs(max(value, key=abs))
        #ax.axis([1.05, -0.05, min(value) - 0.1*max_val, max(value) + 0.1*max_val])



    def plot_me(self):
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1)
        self.my_subplot(ax1, self.voltage, 'voltage / V')
        self.my_subplot(ax2, self.current, 'current / A')
        self.my_subplot(ax3, self.SOC, 'SOC / -')
        self.my_subplot(ax4, self.heat, 'heat / W')
        fig.suptitle("battery output")

        plt.show()

import numpy as np
from fenics import *
from battery_properties import *


class vector_field:

    def __init__(self, eti):
        [num_x, num_y, length_x, length_y] = eti.get_mesh_info()

        self.num_x = int(num_x)
        self.num_y = int(num_y)

        self.length_x = length_x
        self.length_y = length_y


        self.data = np.zeros([self.num_x, num_y, 2], dtype = float)



    def from_fenics(self, fenics_field):

        dx = self.length_x / self.num_x
        dy = self.length_y / self.num_y
        
        for i in range(self.num_x):
            for j in range(self.num_y):
                x = i * dx - self.length_x/2
                y = j * dy - self.length_y/2

                [self.data[i, j, 0], self.data[i, j, 1], tmp] = fenics_field(x, y, 0)
        


class scalar_field:


    def from_eti(self, eti):
        [num_x, num_y, length_x, length_y] = eti.get_mesh_info()
        self.num_x = num_x
        self.num_y = num_y
        self.data = np.zeros(shape=(num_x, num_y))

    def from_vector_field(self, vector_field):

        self.num_x = vector_field.num_x
        self.num_y = vector_field.num_y

        self.data = np.zeros([self.num_x, self.num_y], dtype = float)
        for i in range(self.num_x):
            for j in range(self.num_y):
                self.data[i, j] = sqrt(vector_field.data[i,j,0]**2 + vector_field.data[i,j,1]**2)


    def from_scalar_field(self, other):

        self.num_x = other.num_x
        self.num_y = other.num_y
        self.data = np.zeros([self.num_x, self.num_y], dtype = float)


    def mirror(self):
        new_self = scalar_field()
        new_self.from_scalar_field(self)
        for i in range(self.num_x):
            for j in range(self.num_y):
                jj = self.num_y - j - 1
                new_self.data[i, j] = self.data[i, jj]
        
        return new_self
    

    def __add__(self, other):
        new_self = scalar_field()
        new_self.from_scalar_field(self)
        new_self.data = self.data + other.data

        return new_self
    
    def __mul__(self, scalar):
        new_self = scalar_field()
        new_self.from_scalar_field(self)
        for i in range(self.num_x):
            for j in range(self.num_y):
                new_self.data[i, j] = self.data[i, j] * scalar
        
        return new_self 
    

    def scalar_cell_current(self):
        mirrored_self = self.mirror()
        
        return self + mirrored_self
    
    def sum(self):
        sum = 0
        for i in range(self.num_x):
            for j in range(self.num_y):
                sum += self.data[i, j]
        return sum 
    




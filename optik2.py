# -*- coding: utf-8 -*-
import numpy as np

def g_function(ni,ki,n_back=1.0,k_back=0):
    ans  = (n_back**2 -ni**2 + k_back**2 -ki**2) / (
        (n_back + ni)**2 + (k_back + ki)**2
        )
    return ans
def h_function(ni,ki,n_back=1.0,k_back=0):
    ans  = 2*(n_back*ki - ni*k_back) / ((n_back + ni)**2 + (k_back + ki)**2)
    return ans
def alpha_function(d, k, lambd):
    ans = (2 * np.pi * k * d) / lambd
    return ans
def gamma_function(d, n, lambd):
    ans = (2 * np.pi * n * d) / lambd
    return ans
def read_nk_file():
    file_name = input('File Name:')
    data = np.loadtxt(file_name,dtype=float, skiprows=1) 
    data = data.T
    if len(data) != 3:
        print('3 cols with Wavelength(nm) n k. Try again!')
        return False
    else:
        wl, n, k = data
        step = 0.5
        new_wl = np.arange(280.0,wl.max(),step,dtype=float)
        new_n = np.interp(new_wl,wl,n)
        new_k = np.interp(new_wl,wl,k)
        new_data = np.array([
                list(new_wl),
                list(new_n),
                list(new_k)
            ],dtype=float)
        return new_data


class lego:
    def __init__(self, n, k, name, new=True):
        self.name = name
        if new:
            data = read_nk_file()
            if data:
                self.k, self.n, self.lambd = data[0], data[1], data[2]
            else:
                raise ValueError("Can't upload data")
        else:
            raise "Buscar name en la base, no sé cómo :v"

    def use(self, t, first=False, bottom=False):
        self.thickness = t
        self.alpha = alpha_function(self.thickness, self.k, self.lambd)
        self.gamma = gamma_function(self.thickness, self.n, self.lambd)
        if first:
            self.first = first
            self.bottom = bottom
            self.g = g_function(self.n, self.k)
            self.h = h_function(self.n, self.k)
        if bottom:
            self.first = first
            self.bottom = bottom
            self.g = None
            self.h = None
        return self


class lego_tower():
    def __init__(self, *args):
        self.n_layers = len(args) if len(args) >= 2 else "ERROR"
        if self.n_layers == "ERROR":
            raise "The cell must 2 layers or more"
        self.layers = []
        first = True
        for layer in args:
            thickness = float(input('Thicknes of {}:'.format(layer.name)))
            if first:
                layer.use(thickness, first=first)
                first = False
            else:
                layer.use(thickness, bottom=layer==args[-1])
            self.layers.append(layer)
    def prepare_legos(self):
        for layer in self.layers:
            if layer.first:
                up = layer
                pass
            layer.g = g_function(layer.n, layer.k, up.n, up.k)
            layer.h = h_function(layer.n, layer.k, up.n, up.k)
            up = layer
        return False

    def calc_RT(self):
        self.prepare_legos()
        return False

"""
n = int(input('Number of layers: '))
data = list()
for i in range(n):
    data.append(np.genfromtxt(input('Name of file layer '+str(i+1)+':')))
    if i == 0:
        lambd = data[0][:,][:, 0]
    thickness_layers = float(input('Thickness of file layer '+str(i+1)+': '))
    # thickness_layers = list(float(x) for x in input().split(',')) #for read a LIST INPUT
"""
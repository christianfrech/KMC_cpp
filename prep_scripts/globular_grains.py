import numpy as np 
import math


grain_diam = 10 #angstroms
a = 0.5 #angstroms
c = 1.633 * a #angstroms
a_grain = grain_diam * a
c_grain = grain_diam * c
grain_n_ucells = np.floor(grain_diam / a)
x_dim = 100
y_dim = 100
z_dim = 200

y_spacing = np.sqrt(np.power(grain_n_ucells, 2) - np.power(grain_n_ucells / 2, 2))

a1 = np.array([a_grain/2, -math.sqrt(3) * a_grain/2, 0])
a2 = np.array([a_grain/2, math.sqrt(3) * a_grain/2, 0])
a3 = np.array([0, 0, c_grain])

vecs_along_directions = [5,5,5]
grain_locs = np.zeros((vecs_along_directions[0]*vecs_along_directions[1]*vecs_along_directions[2], 3))



for i in range(vecs_along_directions[0]):
    for j in range(vecs_along_directions[1]):
        for k in range(vecs_along_directions[2]):
            print(f"\nindex: {i*vecs_along_directions[0]*vecs_along_directions[1] + j*vecs_along_directions[0] + k}")
            print(a1 * i + a2 * j + a3 * k)
            grain_locs[i*vecs_along_directions[0]*vecs_along_directions[1] + j*vecs_along_directions[0] + k] = a1 * i + a2 * j + a3 * k

print(f"grain_locs: {grain_locs}")
min = np.min(grain_locs)
grain_locs = grain_locs - np.array([0,1,0]) * min

print(f"grain_locs: {grain_locs}")

'''
### making layer 1 of hcp ###
### row 1 ###
x1_l1 = np.arange(0, x_dim + 1, grain_n_ucells)
y1_l1 = np.arange(0, y_dim + 1, 2 * y_spacing)
z1_l1 = np.arange(0, z_dim + 1,  2 * y_spacing)

### row 2 ###
x1_l1 = np.arange((grain_n_ucells / 2), (xdim - grain_n_ucells / 2 + 1), grain_n_ucells)
y1_l1 = np.arange(y_spacing, y_dim + 1, y_spacing * 2)
z1_l1 = np.arange(0, z_dim + 1, 2 * y_spacing)

### making layer 2 of hcp ###

x1_l2 = np.arange((grain_n_ucells), x_dim + 1, grain_n_ucells)
y1_l2 = np.arange(y_spacing, y_dim + 1, 2 * y_spacing)
z1_l2 = np.arange(y_spacing, z_dim + 1,  2 * y_spacing)

### row 2 ###
x1_l2 = np.arange(0, (xdim + 1), grain_n_ucells)
y1_l2 = np.arange(grain_n_ucells / 2 + y_spacing, y_dim + 1, y_spacing * 2)
z1_l2 = np.arange(y_spacing, z_dim + 1,  2 * y_spacing)
'''

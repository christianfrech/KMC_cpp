import numpy as np 


grain_diam = 10 #angstroms
a = 3.502 #angstroms
grain_n_ucells = np.floor(grain_diam / a)
x_dim = 100
y_dim = 100
z_dim = 200

y_spacing = np.sqrt(np.pow(grain_n_ucells, 2) - np.pow(grain_n_ucells / 2, 2))

a1 = [a/2, -math.sqrt(3)/2, 0]
a2 = [a/2, math.sqrt(3)/2, 0]
a3 = [0, 0, c]


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

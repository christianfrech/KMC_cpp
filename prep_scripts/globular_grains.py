import numpy as np 
import math
import matplotlib.pyplot as plt
import os
import random

def sphere(shape, radius, position):
    """Generate an n-dimensional spherical mask."""
    # assume shape and position have the same length and contain ints
    # the units are pixels / voxels (px for short)
    # radius is a int or float in px
    assert len(position) == len(shape)
    n = len(shape)
    semisizes = (radius,) * len(shape)

    # genereate the grid for the support points
    # centered at the position indicated by position
    grid = [slice(-x0, dim - x0) for x0, dim in zip(position, shape)]
    position = np.ogrid[grid]
    # calculate the distance of all points from `position` center
    # scaled by the radius
    arr = np.zeros(shape, dtype=float)
    for x_i, semisize in zip(position, semisizes):
        # this can be generalized for exponent != 2
        # in which case `(x_i / semisize)`
        # would become `np.abs(x_i / semisize)`
        arr += (x_i / semisize) ** 2

    # the inner part of the sphere will have distance below or equal to 1
    return arr <= 1.0


def idxtoxyz(idxs, x_dim, y_dim, z_dim):
    xyz_idx = np.zeros((3, len(idxs)))
    print(len(idxs))
    xyz_idx[0][:] = idxs[:] % x_dim
    xyz_idx[1][:] = (idxs[:] // (x_dim)) % y_dim
    xyz_idx[2][:] = idxs[:] // (x_dim * y_dim)

    return np.transpose(xyz_idx).astype(int)


def crop_sphere(arr, dim, cell_dims):
    print(f"arr.shape: {arr.shape}")
    print(f"dim {dim}")
    print(f"dim {dim}")
    
    if (dim[0][0] < 0): 
        print("case 1")
        arr = np.delete(arr, np.s_[:(-dim[0][0])], 0)
        print(f"arr: {arr.shape}")
        dim[0][0] = 0

    if (dim[0][1] > cell_dims[0]): 
        print("case 2")
        if (cell_dims[0]-dim[0][1] == 0): slice = -1
        else: slice = cell_dims[0]-dim[0][1]
        arr = np.delete(arr, np.s_[slice:], 0)
        print(f"arr: {arr.shape}")
        dim[0][1] = cell_dims[0]

    if (dim[1][0] < 0): 
        print("case 3")
        arr = np.delete(arr, np.s_[:(-dim[1][0])], 1)
        print(f"arr: {arr.shape}")
        dim[1][0] = 0

    if (dim[1][1] > cell_dims[1]):
        print("case 4")
        if (cell_dims[1]-dim[1][1] == 0): slice = -1
        else: slice = cell_dims[1]-dim[1][1]
        arr = np.delete(arr, np.s_[slice:], 1)
        print(f"arr: {arr.shape}")
        dim[1][1] = cell_dims[1]

    if (dim[2][0] < 0): 
        print("case 5")
        arr = np.delete(arr, np.s_[:(-dim[2][0])], 2)
        print(f"arr: {arr.shape}")
        dim[2][0] = 0

    if (dim[2][1] > cell_dims[2]):
        print("case 6")
        if (cell_dims[2]-dim[2][1] == 0): slice = -1
        else: slice = cell_dims[2]-dim[2][1]
        arr = np.delete(arr, np.s_[slice:], 2)
        print(f"arr: {arr.shape}")
        dim[2][1] = cell_dims[2]

    print(f"arr.shape: {arr.shape}")
    return arr, dim


grain_diam = 25 #angstroms
a = 1 #angstroms
c = 1.633 * a #angstroms
a_grain = grain_diam * a
c_grain = grain_diam * c
grain_n_ucells = np.floor(grain_diam / a)
x_dim = 50
y_dim = 50
z_dim = 50
vec_spacing = 1
vecs_along_directions_start = [-4,-4,-4]
vecs_along_directions_end = [20,20,20]
atomtypes = ["0:vacancy", "1:lithium"]

y_spacing = np.sqrt(np.power(grain_n_ucells, 2) - np.power(grain_n_ucells / 2, 2))

a1 = np.array([a_grain/2, -math.sqrt(3) * a_grain/2, 0])
a2 = np.array([a_grain/2, math.sqrt(3) * a_grain/2, 0])
a3 = np.array([0, 0, c_grain])

shift_layer2 = 1/3*a1 + 2/3*a2 + 1/2*a3 #np.array([])


### drawing centers of spheres in hcp patter ###
num_vecs = (np.floor((np.array(vecs_along_directions_end) - np.array(vecs_along_directions_start)) * np.ones_like(vecs_along_directions_start)/vec_spacing)).astype(int)
print(f"num_vecs: {num_vecs}")
grain_locs = np.zeros(((2*int(num_vecs[0])*int(num_vecs[1])*int(num_vecs[2])), 3), dtype = int)
idx = 0

for i in np.arange(vecs_along_directions_start[0], vecs_along_directions_end[0], vec_spacing):
    for j in np.arange(vecs_along_directions_start[1], vecs_along_directions_end[1], vec_spacing):
        for k in np.arange(vecs_along_directions_start[2], vecs_along_directions_end[2], vec_spacing):
            #print(f"\nidx: {idx}")
            #print(a1 * i + a2 * j + a3 * k)
            grain_locs[idx] = np.floor(a1 * i + a2 * j + a3 * k).astype(int)
            idx += 1
            grain_locs[idx] = np.floor(a1 * i + a2 * j + a3 * k + shift_layer2).astype(int)
            idx += 1

### masking for only points (centers of spheres) which lay inside of the 
### simulation cell ###
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
all_points = np.transpose(np.array(grain_locs), (1,0))

mask = ((all_points[0] <= x_dim) & (all_points[0] >= 0) &
    (all_points[1] <= y_dim) & (all_points[1] >= 0) & 
    (all_points[2] <= z_dim) & (all_points[2] >= 0))

print(f"all_points.shape: {all_points.shape}")
all_points_masked = []
for dimension in all_points: all_points_masked.append(dimension[mask])


### opening output file for writing and writing header info ###
all_points_masked = np.transpose(all_points_masked, (1,0))
file = open("spherical_grains_test.txt", "w")

dim_list = []
for i in range(len(all_points_masked)):
    dim_list.append([int(grain_n_ucells), tuple(all_points_masked[i])])

region_types = [] 
region_coeff = [] 


### generating and assigning random vacancies ###
print("generating and assigning random vacancies")

print("writing lattice to output file")
file.write(f"lattice_dims: {x_dim} {y_dim} {z_dim}\n")
file.write(f"atomtypes: {atomtypes[0]} {atomtypes[1]}\n")
file.write("regions begin\n")

for i in range(len(region_types)):
    if (region_types[i] == "BLOCK"):
        file.write(f"{i+1}: {region_types[i]} xmin:{region_coeff[i][0][0]} xmax:{region_coeff[i][1][0]} ymin:{region_coeff[i][0][1]} ymax:{region_coeff[i][1][1]} zmin:{region_coeff[i][0][2]} zmax:{region_coeff[i][1][2]}\n")

file.write("regions end\n")

### drawing hollow spheres as 1's in boolean array ###
sphere_i = 0
ver_random_vac_xyz = np.ones((x_dim, y_dim, z_dim))
bc_random_vac_xyz = np.ones((x_dim, y_dim, z_dim))

for sphere_dim in dim_list:
    print(f"sphere_i: {sphere_i}")
    print(f"sphere_dim: {sphere_dim}")
    radius = int(sphere_dim[0] / 2)
    #center_radius = sphere_dim[0] - 1

    sphere_arr = sphere((2*radius+1,2*radius+1,2*radius+1), radius, (radius,radius,radius))
    #center_sphere_arr = sphere((2*radius+1,2*radius+1,2*radius+1), 
    #                        center_radius, (center_radius,center_radius,center_radius))


    x_start = int(np.floor(sphere_dim[1][0] - sphere_arr.shape[0]/2))
    x_end = int(np.floor(sphere_dim[1][0] + sphere_arr.shape[0]/2))
    y_start = int(np.floor(sphere_dim[1][1] - sphere_arr.shape[1]/2))
    y_end = int(np.floor(sphere_dim[1][1] + sphere_arr.shape[1]/2))
    z_start = int(np.floor(sphere_dim[1][2] - sphere_arr.shape[2]/2))
    z_end = int(np.floor(sphere_dim[1][2] + sphere_arr.shape[2]/2))
    dims = [[x_start,x_end],[y_start,y_end],[z_start,z_end]]

    ### cutting out center of sphere, leaving outline of sphere ###
    #sphere_arr = (sphere_arr & np.logical_not(center_sphere_arr))
    print(f"hollow_sphere_arr.shape: {sphere_arr.shape}")
    print(f"x: {dims[0][0]}:{dims[0][1]}")
    print(f"y: {dims[1][0]}:{dims[1][1]}") 
    print(f"z: {dims[2][0]}:{dims[2][1]}\n")

    ### eliminating portions of spheres outside of simulation cell ###
    if ((x_start < 0) or (x_end >= x_dim) or 
        (y_start < 0) or (y_end >= y_dim) or 
        (z_start < 0) or (z_end >= z_dim)):
        sphere_arr, dims = crop_sphere(sphere_arr, dims, [x_dim, y_dim, z_dim])
        
    
    #print(np.nonzero(np.logical_not(ver_random_vac_xyz)))
    ver_random_vac_xyz[dims[0][0]:dims[0][1],dims[1][0]:dims[1][1],dims[2][0]:dims[2][1]] = np.logical_not(sphere_arr)
    bc_random_vac_xyz[dims[0][0]:dims[0][1],dims[1][0]:dims[1][1],dims[2][0]:dims[2][1]] = np.logical_not(sphere_arr)

    onegrain_ver_random_vac_xyz = np.ones((x_dim, y_dim, z_dim))
    onegrain_bc_random_vac_xyz = np.ones((x_dim, y_dim, z_dim))
    onegrain_ver_random_vac_xyz[dims[0][0]:dims[0][1],dims[1][0]:dims[1][1],dims[2][0]:dims[2][1]] = np.logical_not(sphere_arr)
    onegrain_bc_random_vac_xyz[dims[0][0]:dims[0][1],dims[1][0]:dims[1][1],dims[2][0]:dims[2][1]] = np.logical_not(sphere_arr)

    ver_nonz_not = np.nonzero(np.logical_not(onegrain_ver_random_vac_xyz))
    ax.scatter(ver_nonz_not[0], ver_nonz_not[1], ver_nonz_not[2], marker='s')



    ### converting array indices and values to strings for file output ###
    #print(np.nonzero(np.logical_not(ver_random_vac_xyz)))
    for x in range(len(ver_random_vac_xyz)):
        for y in range(len(ver_random_vac_xyz[0])):
            for z in range(len(ver_random_vac_xyz[0][0])):
                if (not ver_random_vac_xyz[x][y][z] ):
                    output = f"v {x} {y} {z} {ver_random_vac_xyz[x][y][z]} \n"
                    file.write(output)
                if (not bc_random_vac_xyz[x][y][z]):
                    output = f"bc {x+0.5} {y+0.5} {z+0.5} {bc_random_vac_xyz[x][y][z]} \n"
                    file.write(output)

    sphere_i += 1

file.close()


ax.set_xlim(0,x_dim)
ax.set_ylim(0,y_dim)
ax.set_zlim(0,z_dim)
plt.show()
ax.clear()

ver_nonz = np.nonzero(ver_random_vac_xyz)
ax.scatter(ver_nonz[0], ver_nonz[1], ver_nonz[2], marker='s', alpha=0.025)

ax.set_xlim(0,x_dim)
ax.set_ylim(0,y_dim)
ax.set_zlim(0,z_dim)
plt.show()
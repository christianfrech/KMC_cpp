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


def embed_hexagon_subroutine(sim_cell, locs, size, dims):
    for idx in range(len(locs)):
        grain_loc = locs[idx]
        for k in range(size): # length of grain 
            i3 = k + int(size/2) + grain_loc[2]
            if (i3 >= 0) and (i3 < dims[2]):
                arr = draw_hexagon((size, size, size), int(size/2))
                
                if (k < np.ceil(size/4)): 
                    width = int(4 * k)
                    #print(f"k: {k}, width: {width}")
                    arr = draw_hexagon((size, size), width)
                
                elif (k > np.floor(3*size/4)): 
                    width = int(size - (k - np.floor(3*size/4)) * 4)
                    #print(f"k: {k}, width: {width}")
                    arr = draw_hexagon((size, size), width)
                
                else:
                    width = size 
                    arr = draw_hexagon((size, size), width)
                
                print(f"k: {k}")
                print(f"dims[0]/2: {size/2} dims[1]/2: {size/2} width/2: {arr.shape[0]/2}")
                print(f"shift i: {size/2 - arr.shape[0]/2} shift j: {size/2 - arr.shape[1]/2}")
                print(f"arr.shape: {arr.shape}\n")
                for i in range(len(arr)):
                    for j in range(len(arr[0])):
                        i1 = int(grain_loc[0] + (size/2 - arr.shape[0]/2) + i)
                        i2 = int(grain_loc[1] + (size/2 - arr.shape[1]/2) + j)

                        if (i1 < 0) or (i1 >= dims[0]): pass
                        elif (i2 < 0) or (i2 >= dims[1]): pass
                        elif (arr[i][j] == 1): 
                            #print("check")
                            sim_cell[i1][i2][i3] = 1

    return sim_cell 


### function for generating a hexagon in an arbitrarily-sized box given the sidelength 
### must have (box width >= sidelength*2) and must use odd sidelength*sqrt(3) for odd boxsize
def draw_hexagon(boxsize, sidelen):
    array = np.zeros((2*sidelen, 2*sidelen))
    array_out = np.zeros((int(boxsize[0]), int(boxsize[1])))#int(np.ceil(np.sqrt(2)*sidelen))))
    #print(f"array shape: {array.shape}")
    
    for i in np.arange(int(np.ceil(sidelen/2)), int(np.floor(3*sidelen/2))):
        start_idx = 0
        end_idx = int(np.ceil(np.sqrt(3)*sidelen))
        for j in np.arange(start_idx, end_idx):
            array[i][j] = 1 
    
    for i in np.arange(0, int(np.ceil(sidelen/2))):
        idx = i 
        if (int(np.floor(np.sqrt(3)*sidelen)) % 2 == 0): end_idx = int(np.floor(np.sqrt(3)*sidelen/2)) + 2*idx  - 1
        else: end_idx = int(np.floor(np.sqrt(3)*sidelen/2)) + 2*idx
        start_idx = int(np.ceil(np.sqrt(3)*sidelen/2)) - 2*idx + 1
        
        for j in np.arange(start_idx, end_idx):
            array[i][j] = 1 
    
    for i in np.arange(int(np.floor(3*sidelen/2)), int(2*sidelen)):
        idx = i - int(np.floor(3*sidelen/2))
        start_idx = 2*idx + 2
        end_idx = int(np.floor(np.sqrt(3)*sidelen)) - 2*idx -1
        for j in np.arange(start_idx, end_idx):
            array[i][j] = 1
    
    int_shape = array.shape
    out_shape = array_out.shape
    margin_x = int((out_shape[0]- int_shape[0])/2)
    margin_y = int((out_shape[1]- int_shape[1])/2)
    
    array_out = np.concatenate((np.zeros((array.shape[0], 1)), array), axis=1)
    array_out = np.concatenate((array_out, np.zeros((array.shape[0], 1))), axis=1)
    array_out = np.concatenate((array_out, np.zeros((1, array_out.shape[1]))), axis=0)
    array_out = np.concatenate((np.zeros((1, array_out.shape[1])), array_out), axis=0)
    
    #print(f"array_out.shape: {array_out.shape}")

    return array_out

grain_diam = 10 #angstroms
a = 1 #angstroms
c = 1.633 * a #angstroms
a_grain = grain_diam * a
c_grain = grain_diam * c
grain_n_ucells = np.floor(grain_diam / a)
x_dim = 30
y_dim = 30
z_dim = 30
vec_spacing = 1
vecs_along_directions_start = [-4,-4,-4]
vecs_along_directions_end = [5,5,5]
atomtypes = ["0:vacancy", "1:lithium"]

a1 = np.array([a_grain/2, -math.sqrt(3) * a_grain/2, 0])
a2 = np.array([a_grain/2, math.sqrt(3) * a_grain/2, 0])
a3 = np.array([0, 0, c_grain])

shift_layer2 = 1/3*a1 + 2/3*a2 + 1/2*a3 

### drawing centers of spheres in hcp pattern ###
num_vecs = (np.floor((np.array(vecs_along_directions_end) - np.array(vecs_along_directions_start)) * np.ones_like(vecs_along_directions_start)/vec_spacing)).astype(int)
print(f"num_vecs: {num_vecs}")
grain_locs = np.zeros(((2*int(num_vecs[0])*int(num_vecs[1])*int(num_vecs[2])), 3), dtype = int)
idx = 0

for i in np.arange(vecs_along_directions_start[0], vecs_along_directions_end[0], vec_spacing):
    for j in np.arange(vecs_along_directions_start[1], vecs_along_directions_end[1], vec_spacing):
        for k in np.arange(vecs_along_directions_start[2], vecs_along_directions_end[2], vec_spacing):
            grain_locs[idx] = np.floor(a1 * i + a2 * j + a3 * k).astype(int)
            idx += 1
            grain_locs[idx] = np.floor(a1 * i + a2 * j + a3 * k + shift_layer2).astype(int)
            idx += 1

### masking for only points (centers of spheres) which lay inside of the 
### simulation cell ###
all_points = np.transpose(np.array(grain_locs), (1,0))
mask = ((all_points[0] <= x_dim) & (all_points[0] >= 0) &
    (all_points[1] <= y_dim) & (all_points[1] >= 0) & 
    (all_points[2] <= z_dim) & (all_points[2] >= 0))

print(f"all_points.shape: {all_points.shape}")
all_points_masked = []
for dimension in all_points: all_points_masked.append(dimension[mask])

all_points_masked = np.transpose(all_points_masked, (1,0))
fig = plt.figure()
print(f"all_points_masked.shape: {all_points_masked.shape}")
grain_points = np.transpose(grain_locs, (1,0))
all_points_masked = [[int(np.floor(x_dim/2)),int(np.floor(y_dim/2)),int(np.floor(z_dim/2))]]

ax = fig.add_subplot(111, projection='3d')
sim_cell = np.zeros((x_dim,y_dim, z_dim))
print(f"grain_locs: {grain_locs}")
sim_cell = embed_hexagon_subroutine(sim_cell, all_points_masked, grain_diam, [x_dim, y_dim, z_dim])

arr_nonz = np.nonzero(sim_cell)
ax.scatter(arr_nonz[0], arr_nonz[1], arr_nonz[2])
plt.show()
ax.clear()
'''
ax = fig.add_subplot(111, projection='3d')
arr = generate_grain(10, 8)
arr_nonz = np.nonzero(arr)
ax.scatter(arr_nonz[0], arr_nonz[1], arr_nonz[2])
plt.show()
'''
exit()

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
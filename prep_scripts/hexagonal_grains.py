import numpy as np 
import math
import matplotlib.pyplot as plt

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
        print(f"grain_loc: {grain_loc}")
        for k in range(size): # length of grain 
            i3 = k + grain_loc[2]
            if (i3 >= 0) and (i3 < dims[2]):
                
                if (k < np.ceil(size/4)): 
                    width = int(4 * k)
                    #print(f"k: {k}, width: {width}")
                    arr = draw_hexagon((size, size), int(width/2))
                
                elif (k > np.floor(3*size/4)): 
                    width = int(size - (k - np.floor(3*size/4)) * 4)
                    #print(f"k: {k}, width: {width}")
                    arr = draw_hexagon((size, size), int(width/2))
                
                else:
                    width = size 
                    arr = draw_hexagon((size, size), int(width/2))
                
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
'''
def draw_hexagon_orig(boxsize, sidelen):
    array = np.zeros((2*sidelen, 2*sidelen))
    array_out = np.zeros((int(boxsize[0]), int(boxsize[1]))) # int(np.ceil(np.sqrt(2)*sidelen))))
    
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

    array_out = np.concatenate((np.zeros((array.shape[0], 1)), array), axis=1)
    array_out = np.concatenate((array_out, np.zeros((array.shape[0], 1))), axis=1)
    array_out = np.concatenate((array_out, np.zeros((1, array_out.shape[1]))), axis=0)
    array_out = np.concatenate((np.zeros((1, array_out.shape[1])), array_out), axis=0)

    return array_out
'''

def draw_hexagon(boxsize, sidelen):
    array = np.zeros((int(np.floor(3*sidelen/2)), 2*sidelen))
    array_out = np.zeros((int(boxsize[0]), int(boxsize[1]))) # int(np.ceil(np.sqrt(2)*sidelen))))
    
    for i in np.arange(int(np.ceil(sidelen/2)), int(np.floor(sidelen))):
        start_idx = 0
        end_idx = int(np.ceil(np.sqrt(3)*sidelen))

        for j in np.arange(start_idx, end_idx):
            #print(f"i: {i}, j: {j}")
            array[i][j] = 1 
    
    for i in np.arange(0, int(np.ceil(sidelen/2))):
        idx = i 
        if (int(np.floor(np.sqrt(3)*sidelen)) % 2 == 0): end_idx = int(np.floor(np.sqrt(3)*sidelen/2)) + 2*idx  - 1
        else: end_idx = int(np.floor(np.sqrt(3)*sidelen/2)) + 2*idx
        start_idx = int(np.ceil(np.sqrt(3)*sidelen/2)) - 2*idx + 1
        
        for j in np.arange(start_idx, end_idx):
            #print(f"i: {i}, j: {j}")
            array[i][j] = 1 
    
    for i in np.arange(int(sidelen), int(np.floor(3*sidelen/2))):
        idx = i - int(np.floor(sidelen))
        start_idx = 2*idx + 2
        end_idx = int(np.floor(np.sqrt(3)*sidelen)) - 2*idx -1
        
        for j in np.arange(start_idx, end_idx):
            #print(f"i: {i}, j: {j}")
            array[i][j] = 1

    array_out = np.concatenate((np.zeros((array.shape[0], 1)), array), axis=1)
    array_out = np.concatenate((array_out, np.zeros((array.shape[0], 1))), axis=1)
    array_out = np.concatenate((array_out, np.zeros((1, array_out.shape[1]))), axis=0)
    array_out = np.concatenate((np.zeros((1, array_out.shape[1])), array_out), axis=0)

    return array_out


grain_diam = 30 #angstroms

a = 1 #angstroms
c = 1.633 * a #angstroms
a_grain = grain_diam * a
c_grain = grain_diam * c
grain_n_ucells = np.floor(grain_diam / a)

x_dim = 100
y_dim = 100
z_dim = 100

vecs_along_directions_start = [-4,-4,-4]
vecs_along_directions_end = [20,20,20]
atomtypes = ["0:vacancy", "1:lithium"]

vec_spacing = 1
xscale = 0.75
yscale = 0.8
zscale = 0.9
a1 = np.array([a_grain/2 * xscale, -math.sqrt(3) * a_grain/2 * yscale, 0])
a2 = np.array([a_grain/2 * xscale, math.sqrt(3) * a_grain/2 * yscale, 0])
a3 = np.array([0, 0, c_grain * zscale])

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

mask = ((all_points[0] <= (x_dim + grain_diam)) & (all_points[0] >= -grain_diam) &
    (all_points[1] <= (y_dim + grain_diam)) & (all_points[1] >= -grain_diam) & 
    (all_points[2] <= (z_dim + grain_diam)) & (all_points[2] >= -grain_diam))

all_points_masked = []
for dimension in all_points: all_points_masked.append(dimension[mask])

all_points_masked = np.transpose(all_points_masked, (1,0))
fig = plt.figure()
grain_points = np.transpose(grain_locs, (1,0))

ax = fig.add_subplot(111, projection='3d')
sim_cell = np.zeros((x_dim,y_dim, z_dim))
print(f"grain_locs: {grain_locs}")
sim_cell = embed_hexagon_subroutine(all_points_masked, grain_diam, [x_dim, y_dim, z_dim])

arr_nonz = np.nonzero(np.logical_not(sim_cell))
#arr_nonz = np.nonzero(sim_cell)
ax.scatter(arr_nonz[0], arr_nonz[1], arr_nonz[2], alpha = 0.1)
print(f"occupancy: {len(arr_nonz[0])/(x_dim*y_dim*z_dim)}")
plt.show()
ax.clear()

### opening output file for writing and writing header info ###
file = open("hexagonal_grains_test.txt", "w")

### converting array indices and values to strings for file output ###
#print(np.nonzero(np.logical_not(ver_random_vac_xyz)))
arr_nonz = np.transpose(arr_nonz, (1,0))
print(f"occupancy: {len(arr_nonz)/(x_dim *y_dim *z_dim)}")
print(arr_nonz.shape)
for nonz in arr_nonz:
    output = f"v {int(nonz[0])} {int(nonz[1])} {int(nonz[2])} \n"
    file.write(output)
    output = f"bc {int(nonz[0])} {int(nonz[1])} {int(nonz[2])} \n"
    file.write(output)
    
    '''
    output = f"v {x} {y} {z} \n" # {sim_cell[x][y][z]} \n"
    file.write(output)
    output = f"bc {x+0.5} {y+0.5} {z+0.5} \n" #  {sim_cell[x][y][z]} \n"
    file.write(output)
    '''
 
file.close()
import os
import numpy as np
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



file = open("test_sphere.txt", "w")
x_dim = 130
y_dim = 130
z_dim = 130
num_of_rand = 2
step = 2
cluster_num = 1
geo = "bcc"
atomtypes = ["0:vacancy", "1:lithium"]
region_bounds = []
all_vacancies = np.ones((2, x_dim, y_dim, z_dim))
dim_list = [[5,(55,55,115)]]

region_types = [] 
region_coeff = [] 

def idxtoxyz(idxs, x_dim, y_dim, z_dim):
    xyz_idx = np.zeros((3, len(idxs)))
    print(len(idxs))
    xyz_idx[0][:] = idxs[:] % x_dim
    xyz_idx[1][:] = (idxs[:] // (x_dim)) % y_dim
    xyz_idx[2][:] = idxs[:] // (x_dim * y_dim)

    return np.transpose(xyz_idx).astype(int)

### generating and assigning random vacancies ###
print("generating and assigning random vacancies")

print("writing lattice to output file")
file.write(f"lattice_dims: {x_dim} {y_dim} {z_dim}\n")
file.write(f"atomtypes: {atomtypes[0]} {atomtypes[1]}\n")
file.write(f"num_regions: {len(region_types)}\n")
file.write("regions begin\n")

for i in range(len(region_types)):
    #if (region_types[i] == "GB"):
        #file.write(f"{i+1}: {region_types[i]} x:{region_coeff[i][0]} y:{region_coeff[i][1]} z:{region_coeff[i][2]} xshift:{region_shift[i][0]} yshift:{region_shift[i][1]} zshift:{region_shift[i][2]}\n")
    if (region_types[i] == "BLOCK"):
        file.write(f"{i+1}: {region_types[i]} xmin:{region_coeff[i][0][0]} xmax:{region_coeff[i][1][0]} ymin:{region_coeff[i][0][1]} ymax:{region_coeff[i][1][1]} zmin:{region_coeff[i][0][2]} zmax:{region_coeff[i][1][2]}\n")

file.write("regions end\n")

for sphere_dim in dim_list:
    print(f"sphere_dim: {sphere_dim}")
    radius = sphere_dim[0]

    sphere_arr = sphere((2*radius+1,2*radius+1,2*radius+1), radius, (radius,radius,radius))
    print(sphere_arr.shape)

    x_start = int(sphere_dim[1][0] - sphere_arr.shape[0]/2)
    x_end = int(sphere_dim[1][0] + sphere_arr.shape[0]/2)
    y_start = int(sphere_dim[1][1] - sphere_arr.shape[1]/2)
    y_end = int(sphere_dim[1][1] + sphere_arr.shape[1]/2)
    z_start = int(sphere_dim[1][2] - sphere_arr.shape[2]/2)
    z_end = int(sphere_dim[1][2] + sphere_arr.shape[2]/2)

    ver_random_vac_xyz = np.ones((x_dim, y_dim, z_dim))

    ver_random_vac_xyz[x_start:x_end,y_start:y_end,z_start:z_end] = np.logical_not(sphere_arr[:][:][:])

    all_vacancies[0,x_start:x_end,y_start:y_end,z_start:z_end] = np.logical_not(sphere_arr[:][:][:])
    all_vacancies[1,x_start:x_end,y_start:y_end,z_start:z_end] = np.logical_not(sphere_arr[:][:][:])
    all_vacancies[1,(x_start-1):(x_end-1),(y_start):(y_end),(z_start):(z_end)] = np.logical_not(sphere_arr[:][:][:])
    all_vacancies[1,(x_start-1):(x_end-1),(y_start-1):(y_end-1),(z_start):(z_end)] = np.logical_not(sphere_arr[:][:][:])
    all_vacancies[1,(x_start-1):(x_end-1),(y_start-1):(y_end-1),(z_start-1):(z_end-1)] = np.logical_not(sphere_arr[:][:][:])
    all_vacancies[1,(x_start-1):(x_end-1),(y_start):(y_end),(z_start-1):(z_end-1)] = np.logical_not(sphere_arr[:][:][:])
    all_vacancies[1,(x_start):(x_end),(y_start-1):(y_end-1),(z_start-1):(z_end-1)] = np.logical_not(sphere_arr[:][:][:])
    all_vacancies[1,(x_start):(x_end),(y_start-1):(y_end-1),(z_start):(z_end)] = np.logical_not(sphere_arr[:][:][:])
    all_vacancies[1,(x_start):(x_end),(y_start):(y_end),(z_start-1):(z_end-1)] = np.logical_not(sphere_arr[:][:][:])


    ### converting array indices and values to strings for file output ###
    print(np.nonzero(np.logical_not(ver_random_vac_xyz)))
    for x in range(len(ver_random_vac_xyz)):
        for y in range(len(ver_random_vac_xyz[0])):
            for z in range(len(ver_random_vac_xyz[0][0])):
                if (not all_vacancies[0][x][y][z] ):
                    output = f"v {x} {y} {z} {all_vacancies[0][x][y][z]} \n"#{cluster_num} \n"
                    file.write(output)

                if (not all_vacancies[1][x][y][z] ):
                    output = f"bc {x + 0.5} {y + 0.5} {z + 0.5} {all_vacancies[1][x][y][z]} \n"# {cluster_num} \n"
                    file.write(output)

    cluster_num += 1


file.close()

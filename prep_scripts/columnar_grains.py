import numpy as np
import random
import matplotlib.pyplot as plt
#import pyvoro
import math
from scipy.spatial import Voronoi, voronoi_plot_2d
import shapely



def FourD_idxs(start_vec, end_vec):

    idxs = []
    print("running")
    print(f"start_vec: {start_vec}")
    print(f"end_vec: {end_vec}")
    for i in np.arange(start_vec[0], end_vec[0]+1):
        for j in np.arange(start_vec[1], end_vec[1]+1):
            for k in np.arange(start_vec[2], end_vec[2]+1):
                for l in np.arange(start_vec[3], end_vec[3]+1):
                    print([i, j, k, l]) 
                    idx = [i, j, k, l]
                    idxs.append(idx)

    return np.array(idxs)

def trapez(y,y0,w):
    return np.clip(np.minimum(y+1+w/2-y0, -y+1+w/2+y0),0,1)

def weighted_line(r0, c0, r1, c1, w, rmin=0, rmax=np.inf):
    # The algorithm below works fine if c1 >= c0 and c1-c0 >= abs(r1-r0).
    # If either of these cases are violated, do some switches.
    if abs(c1-c0) < abs(r1-r0):
        # Switch x and y, and switch again when returning.
        xx, yy, val = weighted_line(c0, r0, c1, r1, w, rmin=rmin, rmax=rmax)
        return (yy, xx, val)

    # At this point we know that the distance in columns (x) is greater
    # than that in rows (y). Possibly one more switch if c0 > c1.
    if c0 > c1:
        return weighted_line(r1, c1, r0, c0, w, rmin=rmin, rmax=rmax)

    # The following is now always < 1 in abs
    slope = (r1-r0) / (c1-c0)

    # Adjust weight by the slope
    w *= np.sqrt(1+np.abs(slope)) / 2

    # We write y as a function of x, because the slope is always <= 1
    # (in absolute value)
    x = np.arange(c0, c1+1, dtype=float)
    y = x * slope + (c1*r0-c0*r1) / (c1-c0)

    # Now instead of 2 values for y, we have 2*np.ceil(w/2).
    # All values are 1 except the upmost and bottommost.
    thickness = np.ceil(w/2)
    yy = (np.floor(y).reshape(-1,1) + np.arange(-thickness-1,thickness+2).reshape(1,-1))
    xx = np.repeat(x, yy.shape[1])
    vals = trapez(yy, y.reshape(-1,1), w).flatten()

    yy = yy.flatten()

    # Exclude useless parts and those outside of the interval
    # to avoid parts outside of the picture
    mask = np.logical_and.reduce((yy >= rmin, yy < rmax, vals > 0))

    return (yy[mask].astype(int), xx[mask].astype(int), vals[mask])

def check(list1, val):
    return(all(x >= val for x in list1))

def greater_than(list1, val):
    for elem in list1:
        if elem > val: return True
    return False

def less_than(list1, val):
    for elem in list1:
        if elem < val: return True
    return False


dims = [[0,50],[0,32],[0,64]]
x = 5
y_list = [x+1,x,x+1,x,x+1]
num_grains = [len(y_list),y_list,4]
points = np.zeros((np.sum(num_grains[1]), 2))
point_check = np.zeros((dims[0][1] * dims[1][1], 2))
arr = np.zeros((dims[0][1], dims[1][1], dims[2][1]))
plane_val = set([0,1,2])

for x in range(dims[0][1]):
    for y in range(dims[1][1]):
        for z in range(dims[2][1]):
            point_check[(x*dims[1][1] + y)] = [x,y]

### random distribution of points for tesselation ##
'''
for i in range(num_grains[0]):
    for j in range(num_grains[1]):
        for k in range(num_grains[2]):
            grain_num = i*num_grains[1]*num_grains[2]+j*num_grains[2]+k
            points[grain_num][0] = random.randint(0, dims[0][1]-1)
            points[grain_num][1] = random.randint(0, dims[1][1]-1)
'''

### even distribution of points for tesselation ##
grain_num = 0
print(f"num_grains[1].shape: {np.array(num_grains[1]).shape}")

i_step = dims[0][1]/(num_grains[0]-2)
print(f"i_arr: {np.arange(-i_step/2, dims[0][1] + (3/2*i_step) + 1, i_step)}")

for i in range(num_grains[0]):
    i_val = np.arange(-i_step/2, dims[0][1] + (3/2*i_step) + 1, i_step)[i]
    j_step = dims[1][1]/(num_grains[1][i]-2)
    #print(f"j_arr: {np.arange(-j_step/2, (dims[1][1] + (3/2*j_step)+ 1), j_step)}")
    for j in range(num_grains[1][i]):
        j_val = np.arange(-j_step/2, dims[1][1] + (3/2*j_step) + 1, j_step)[j]
        points[grain_num][0] = i_val
        points[grain_num][1] = j_val
        grain_num += 1

for x in range(dims[0][1]):
    for y in range(dims[1][1]):
        point_check[(x*dims[1][1] + y)] = [x,y]

#print(f"points: {points}")

vor = Voronoi(points)
all_points = np.array([])
i = 0

# for each Voronoi cell, plot all the faces of the corresponding polygon
vertices_list = vor.vertices
print(f"vor.ridge_vertices.shape: {np.array(vor.ridge_vertices).shape}")
print(f"vertices_list.shape: {np.array(vertices_list).shape}")


for ridge_idxs in vor.ridge_vertices:
    # the vertices are the corner points of the Voronoi cell
    # creating parallelapiped to represent GB
    
    if not check(ridge_idxs, 0): continue
    ver_coords = vertices_list[ridge_idxs]
    unique, unique_inverse = np.unique(ver_coords, return_inverse=True, axis=0)
    unique = unique[unique_inverse]
    
    if not np.array_equal(ver_coords, unique): continue

    flat_ver_coords = [item for row in ver_coords for item in row]

    if greater_than(flat_ver_coords, 55) or not check(flat_ver_coords, -5): continue

    tup_out = weighted_line(ver_coords[0][0], ver_coords[0][1], ver_coords[1][0], ver_coords[1][1], 0.25, rmin=0, rmax=np.inf)
    tup_out = np.array([tup_out[0], tup_out[1]])
    line_coords = np.transpose(tup_out, (1,0))
    line_coords = np.array(line_coords)
    z_coords = np.arange(dims[2][0], dims[2][1], 1)
    z_coords = np.tile(z_coords, line_coords.shape[0])
    line_coords = np.repeat(line_coords, dims[2][1], axis = 0)
    z_coords = np.transpose(np.array([z_coords]), (1,0))
    line_coords = np.append(line_coords, z_coords, axis = 1)

    if (len(all_points.shape) == 1): all_points = line_coords
    else: 
        all_points = np.append(all_points, line_coords, axis = 0)

reg_idx = 0
points_in = []

for reg in vor.regions:
    point_idx = 0
    points_in_reg = np.zeros((len(point_check) , 2))
    print(f"reg: {reg}")
    if less_than(reg, 0): continue
    coords = vertices_list[reg]
    if (len(coords) <= 2): continue
    reg_idx += 1
    coords_tuple = tuple([tuple(coord) for coord in coords])
    polygon = shapely.Polygon(coords_tuple)
    
    for check_point in point_check: 
        point = shapely.geometry.Point(check_point[0], check_point[1])
        if point.within(polygon): 
            points_in_reg[point_idx] = check_point
            point_idx += 1

    np.array(points_in).resize((point_idx, 2))
    points_in.append(points_in_reg)

print(f"points_in: {[points_in[0], points_in[1], np.ones_like(points_in[0]) * 32]}")

fig = voronoi_plot_2d(vor)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
all_points = np.transpose(np.array(all_points), (1,0))
ax.scatter(all_points[0], all_points[1], all_points[2], marker='s', alpha=0.025)

partition_color = "orange"
all_points = np.transpose(all_points, (1,0))
for i in range(len(points_in)):
    points = np.transpose(points_in[i], (1,0))
    rand = random.randint(0, 4)
    print(f"i: {i}")
    print(f"rand: {rand}")
    print(f"points.shape: {points.shape}")
    print(f"np.array([points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 2)]).shape: {np.array([points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 2)]).shape}")
    print(f"all_points.shape: {all_points.shape}\n")

    if (rand%4 == 1):
        ax.scatter(points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 2), marker='s', alpha=0.5, color = partition_color)
        plane = np.array([points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 2)])
        all_points = np.append(all_points, np.transpose(plane, (1,0)), axis = 0)
    
    if (rand%4 == 2):
        ax.scatter(points[0], points[1], (2 * np.ones_like(points[0]) * dims[2][1] / 3), marker='s', alpha=0.5, color = partition_color)
        plane = np.array([points[0], points[1], (2 * np.ones_like(points[0]) * dims[2][1] / 3)])
        all_points = np.append(all_points, np.transpose(plane, (1,0)), axis = 0)
        ax.scatter(points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 3), marker='s', alpha=0.5, color = partition_color)
        plane = np.array([points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 3)])
        all_points = np.append(all_points, np.transpose(plane, (1,0)), axis = 0)

    if (rand%4 == 3):
        ax.scatter(points[0], points[1], (3 * np.ones_like(points[0]) * dims[2][1] / 4), marker='s', alpha=0.5, color = partition_color)
        plane = np.array([points[0], points[1], (3 * np.ones_like(points[0]) * dims[2][1] / 4)])
        all_points = np.append(all_points, np.transpose(plane, (1,0)), axis = 0)
        ax.scatter(points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 2), marker='s', alpha=0.5, color = partition_color)
        plane = np.array([points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 2)])
        all_points = np.append(all_points, np.transpose(plane, (1,0)), axis = 0)
        ax.scatter(points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 4), marker='s', alpha=0.5, color = partition_color)
        plane = np.array([points[0], points[1], (np.ones_like(points[0]) * dims[2][1] / 4)])
        all_points = np.append(all_points, np.transpose(plane, (1,0)), axis = 0)

all_points = np.unique(all_points, axis=0)
print(f"all_points.shape: {all_points.shape}")

file = open("gb_sites.txt", "w")
print(f"all_points: {all_points}")
for point in all_points:
    if (point[2]%2 == 0):
        output = f"v {point[0]} {point[1]} {point[2]} \n"
        file.write(output)
    else:
        output = f"bc {point[0]} {point[1]} {point[2]} \n"
        file.write(output)

plt.show()

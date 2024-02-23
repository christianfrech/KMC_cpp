import numpy as np
from scipy.spatial import ConvexHull
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.cm
import pyvoro
import math


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

def in_poly_hull_multi(poly, points):
    hull = ConvexHull(poly)
    res = []

    for p in points:
        new_hull = ConvexHull(np.concatenate((poly, [p])))
        res.append(np.array_equal(new_hull.vertices, hull.vertices))

    return res

dims = [[0,50],[0,32],[0,64]]
num_grains = [2,2,2]
point_dims = num_grains
points = np.zeros((num_grains[0] * num_grains[1] * num_grains[2], 3))
point_check = np.zeros((dims[0][1] * dims[1][1] * dims[2][1], 3))
arr = np.zeros((dims[0][1], dims[1][1], dims[2][1]))
plane_val = set([0,1,2])

for i in range(num_grains[0]):
    for j in range(num_grains[1]):
        for k in range(num_grains[2]):
            grain_num = i*num_grains[1]*num_grains[2]+j*num_grains[2]+k
            points[grain_num][0] = random.randint(0, dims[0][1]-1)
            points[grain_num][1] = random.randint(0, dims[1][1]-1)
            points[grain_num][2] = random.randint(0, dims[2][1]-1)

for x in range(dims[0][1]):
    for y in range(dims[1][1]):
        for z in range(dims[2][1]):
            point_check[(x*dims[1][1]*dims[2][1] + y*dims[2][1] + z)] = [x,y,z]

voronoi = pyvoro.compute_voronoi(points,dims,6)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
rng = np.random.default_rng(11)
all_points = np.array([])
polygons = []
i = 0
j = 0
# for each Voronoi cell, plot all the faces of the corresponding polygon
print(f"voronoi len: {len(voronoi)}")
for vnoicell in voronoi:
    print(f"i: {i}")
    faces = []
    # the vertices are the corner points of the Voronoi cell
    vertices = np.array(vnoicell['vertices'])
    # cycle through all faces of the polygon
    print(f"vnoicell len: {len(vnoicell['faces'])}")

    for face in vnoicell['faces']:
        print(f"j: {j}")
        totals = {}
        faces.append(vertices[np.array(face['vertices'])])
        # creating parallelapiped to represent GB
        vertices_transp = np.floor(np.transpose(vertices[np.array(face['vertices'])], (1,0)))
        ones = 2*np.ones_like(vertices_transp)
        print(f"shape: {vertices_transp.shape}")
        print(f"{vertices_transp}") 
        print(f"xset: {set(vertices_transp[0][:])}, yset: {set(vertices_transp[:][1])}, zset: {set(vertices_transp[:][2])}")

        # eliminating polygons which are actually lines (collinear four points)
        if (((len(set(vertices_transp[0][:])) == 1) and (len(set(vertices_transp[1][:])) == 1)) or 
            ((len(set(vertices_transp[0][:])) == 1) and (len(set(vertices_transp[2][:])) == 1)) or 
            ((len(set(vertices_transp[1][:])) == 1) and (len(set(vertices_transp[2][:])) == 1))): continue

        # checking that there are multiple coordinates along each dim
        # and creating prism accordingly
        if (len(set(vertices_transp[2][:])) == 1): 
            print("case 1")
            ones[:2] = 0
        elif (len(set(vertices_transp[1][:])) == 1):
            print("case 2")
            ones[0] = 0
            ones[2] = 0
        elif (len(set(vertices_transp[0][:])) == 1):
            print("case 3")
            ones[1] = 0
            ones[2] = 0

        print(f"ones: {ones}")

        ceil = vertices_transp + ones
        floor = vertices_transp
        ceil = np.transpose(ceil, (1,0))
        floor = np.transpose(floor, (1,0))
        prism = np.append(ceil, floor, axis=0)
        print(f"prism: {prism}")
        # eliminating polygons which are lines due to replica points
        res = list(set(map(lambda i: tuple(sorted(i)), prism)))
        # checking to see which points in 3D mesh lie in GB, appending to list  
        res = in_poly_hull_multi(prism, point_check)

        print(f"point_check[res]: {point_check[res]}\n")
        if (point_check[res].size == 0): continue

        if (plane_val == set(point_check[res][0])): continue
        if (plane_val == set(point_check[res][1])): continue
        if (plane_val == set(point_check[res][2])): continue

        if (len(all_points.shape) == 1): all_points = point_check[res]
        else: 
            all_points = np.append(all_points, point_check[res], axis = 0)
            polygons.append(point_check[res])
            print(f"all_points: {all_points}\n")

        j += 1
    # join the faces into a 3D polygon
    polygon = Poly3DCollection(faces, alpha=0.5, 
                               facecolors=rng.uniform(0,1,3),
                               linewidths=0.5,edgecolors='black')

    ax.add_collection3d(polygon)
    i += 1
    


file = open("gb_sites.txt", "w")

for point in all_points:
    if (point[2]%2 == 0):
        output = f"v {point[0]} {point[1]} {point[2]} \n"
        file.write(output)
    else:
        output = f"bc {point[0]} {point[1]} {point[2]} \n"
        file.write(output)


file.close()
nonzeros = np.nonzero(arr)

ax.set_xlim(dims[0][0], dims[0][1])
ax.set_ylim(dims[1][0], dims[1][1])
ax.set_zlim(dims[2][0], dims[2][1])
plt.savefig("test_viorni.png")
ax.clear()


for polyg in polygons:
    polyg = np.transpose(np.array(polyg), (1,0))
    print(f"polyg: {polyg}")
    if (plane_val == set(polyg[0])): continue
    if (plane_val == set(polyg[1])): continue
    if (plane_val == set(polyg[2])): continue
    
    ax.scatter(polyg[0], polyg[1], polyg[2], marker='s', alpha=0.25)

ax.set_xlim(dims[0][0], dims[0][1])
ax.set_ylim(dims[1][0], dims[1][1])
ax.set_zlim(dims[2][0], dims[2][1])
plt.show()
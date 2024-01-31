from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt 
import numpy as np 
from pylab import *
import copy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import math
import scipy as sci
import os
import sys

np.set_printoptions(threshold=sys.maxsize)

def parse_line(line):
    toks = line.split(" ")
    lattice_pos = (toks[0])
    x = (toks[1])
    y = (toks[2])
    z = (toks[3])
    return lattice_pos,x,y,z


def read_vac_file(filename, size, iterations):
    output = np.zeros((size, 3))
    read_file = open(filename, 'r')
    lines = [line for line in read_file]

    idx = 0
    for line in lines:
        lattice_pos,x,y,z = parse_line(line)
        lattice_pos = int(lattice_pos)
        x = int(x)
        y = int(y)
        z = int(z)
        if (lattice_pos == 0):
            output[idx][0] = (x*2) % 200
            output[idx][1] = (y*2) % 200
            output[idx][2] = (z*2) % 200

        else: 
            output[idx][0] = (x*2 + 1) % 200
            output[idx][1] = (y*2 + 1) % 200
            output[idx][2] = (z*2 + 1) % 200

        idx += 1

    return output

size = 200
iterations = 75
j=0
dots_colors = ["lime", "blue", "red", "orange"]
lines_colors = ["lime", "blue", "red", "orange"]



pwd = os.getcwd()
dir = "counts_output"
Ts = [300, 500, 700, 900]
all_vacancies = np.zeros((1, size*iterations, 3))
rootdir = os.path.join(pwd, dir)
rootdir = os.path.join(rootdir, "vacs")

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        #if (dir == "varied_rates_counts"):
        print(dir)
        print(file)
        filepath = os.path.join(rootdir, file)
        T = int(file.split("_")[4][:-8])
        partition = int(file.split("_")[3])
        idx = int(file.split("_")[2])
        coords = read_vac_file(filepath, size, iterations)
        print(f"idx: {idx}")
        start_idx = (idx) * size
        end_idx = (idx+1) * size
        print(f"start idx: {start_idx}")
        print(f"end idx: {end_idx}")
        print(f"partition: {partition}")
        if (idx < iterations) and (partition == 4) and (T == 0):
            all_vacancies[0, start_idx:end_idx, :] = coords

all_vacancies = np.transpose(all_vacancies, (0,2,1))
all_vacancies = np.mod(all_vacancies + np.ones_like(all_vacancies) * 100, np.ones_like(all_vacancies) * 200)
all_vacancies = np.transpose(all_vacancies, (0,2,1))

fig = plt.figure(figsize=(10, 10)) 
ax = fig.add_subplot(111, projection='3d')

for i in range(len(all_vacancies)):
    print(f"i: {i}")
    
    new_bins_x = np.arange(math.floor(np.min(all_vacancies[i][0])), math.ceil(np.max(all_vacancies[i][0])+1), 1)
    new_bins_y = np.arange(math.floor(np.min(all_vacancies[i][1])), math.ceil(np.max(all_vacancies[i][1])+1), 1)
    new_bins_z = np.arange(math.floor(np.min(all_vacancies[i][2])), math.ceil(np.max(all_vacancies[i][2])+1), 1)
    
    print(all_vacancies.shape)
    print(new_bins_x.shape)
    print(new_bins_y.shape)
    print(new_bins_z.shape)

    counts, bins = np.histogramdd(all_vacancies[i], [new_bins_x, new_bins_y, new_bins_z])
    norm_arr = np.ones_like(counts) * size * iterations
    counts_normed = np.divide(counts, norm_arr) # setting color bar 

    #color_map = cm.ScalarMappable(cmap=cm.plasma)
    #color_map.set_array(counts_normed)
    color_map = counts_normed
    #print(color_map)
    

    # creating figures 
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #all_vacancies = np.transpose(all_vacancies, (0,2,1))
    cmplasma = plt.get_cmap("plasma")
    
    # creating the heatmap
    print(f"counts_normed.shape: {counts_normed.shape}")
    print(f"new_bins_x[:-1].shape: {new_bins_x[:-1].shape}")
    print(f"new_bins_y[:-1].shape: {new_bins_y[:-1].shape}")
    print(f"new_bins_z[:-1].shape: {new_bins_z[:-1].shape}")
    print(f"new_bins_x.shape: {new_bins_x}")
    print(f"new_bins_y.shape: {new_bins_y}")
    print(f"new_bins_z.shape: {new_bins_z}")
    bin_x = np.repeat([new_bins_x[:-1]], len(new_bins_y[:-1]), axis=0)
    bin_x = np.repeat([bin_x], len(new_bins_z[:-1]), axis=0)
    bin_x = np.transpose(bin_x, (2,1,0))
    #print( f"bin_x: {bin_x} " )
    print( f"bin_x.shape: {bin_x.shape} " )
    bin_y = np.repeat([new_bins_y[:-1]], len(new_bins_z[:-1]), axis=0)
    bin_y = np.repeat([bin_y], len(new_bins_x[:-1]), axis=0)
    bin_y = np.transpose(bin_y, (0,2,1))
    #print( f"bin_y: {bin_y} " )
    print( f"bin_y.shape: {bin_y.shape} " )
    bin_z = np.repeat([new_bins_z[:-1]], len(new_bins_y[:-1]), axis=0)
    bin_z = np.repeat([bin_z], len(new_bins_x[:-1]), axis=0)
    #bin_x = np.transpose(bin_x, (2,1,0))
    #print( f"bin_z: {bin_z} " )
    print( f"bin_z.shape: {bin_z.shape} " )

    ax.scatter(bin_x, bin_y, bin_z, counts_normed, marker='s', s=25, alpha=0.4, c=color_map)
    #plt.colorbar(color_map)
    ax.set_xlabel("(z) lattice position")
    ax.set_ylabel("x/x_0")
    plt.savefig(f"Plot3D_bigtest.png")
    plt.cla()

import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
import time
import math
import scipy as sci
import os
import sys
import imageio

np.set_printoptions(threshold=sys.maxsize)

def parse_line(line):
    toks = line.split(" ")
    lattice_pos = (toks[0])
    x = (toks[1])
    y = (toks[2])
    z = (toks[3])
    cluster_num = (toks[4])
    return lattice_pos,x,y,z,cluster_num


def read_cluster_file(filename, size, dims, num_clusters):
    output = np.zeros((size, 3))
    read_file = open(filename, 'r')
    lines = [line for line in read_file]

    start_idxs = np.zeros((num_clusters + 1))

    for line in lines:
        lattice_pos,x,y,z,cluster_num = parse_line(line)
        x = int(float(x))
        y = int(float(y))
        z = int(float(z))
        cluster_num = int(cluster_num)
        start_idxs[(cluster_num)] += 1

    idxs = start_idxs.astype(int)
    
    for line in lines:
        lattice_pos,x,y,z,cluster_num = parse_line(line)
        x = int(float(x))
        y = int(float(y))
        z = int(float(z))
        cluster_num = int(cluster_num)

        if (lattice_pos == "v"):
            clus_idx = idxs[(cluster_num - 1)]
            output[clus_idx][0] = (x*2) % (2*dims[0])
            output[clus_idx][1] = (y*2) % (2*dims[1])
            output[clus_idx][2] = (z*2) % (2*dims[2])

        else: 
            clus_idx = idxs[(cluster_num - 1)]
            output[clus_idx][0] = (x*2 + 1) % (2*dims[0])
            output[clus_idx][1] = (y*2 + 1) % (2*dims[1])
            output[clus_idx][2] = (z*2 + 1) % (2*dims[2])

        idxs[(cluster_num - 1)] += 1

    for i in range(len(idxs)):
        end_idx = idxs[i]
        output[i] = output[i][:end_idx]

    return start_idxs, output


size = 6388
iterations = 1
maxtime = 28100000
start = 100000

times = np.arange(start, maxtime, start, dtype=int)
print(f"times: {times}")
lines_colors = ["blue", "red", "orange", "green"]
dims = [100, 100, 100]
idx=0
fig, ax = plt.subplots()
j=0
num_clusters = 2
dir_names = ["twovoids_4NN"]
data_names = []
all_vacancies = np.zeros((len(dir_names),iterations, len(times), size, 3))
all_times = np.zeros((len(dir_names),iterations,len(times)), dtype = float)
all_start_idxs = np.zeros((len(dir_names),iterations,len(times),(num_clusters+1)))


rootdir = os.getcwd()


for dir in filter(os.path.isdir, os.listdir(os.getcwd())):
    if (dir in dir_names):
        mid_filepath = os.path.join(rootdir, dir)
        mid_filepath = os.path.join(mid_filepath, "cluster_nums")
        data_names.append(dir)
        print(mid_filepath)

        for subdir, dirs, files in os.walk(mid_filepath):
            for file in files:
                    
                filepath = os.path.join(mid_filepath, file)
                line_split = file.split("_")
                iteration = int(line_split[2])
                step = int(line_split[3])
                step_time = float(line_split[4])
                
                if (step_time <= maxtime):
                    idx = np.where(times == step)[0]
                    if (len(idx) == 0):
                        continue
                    
                    start_idxs, coords = read_cluster_file(filepath, size, dims, num_clusters)
                    all_vacancies[j][iteration][idx] = coords
                    all_times[j][iteration][idx] = step_time
                    all_start_idxs[j][iteration][idx] = start_idxs
                    idx += 1

        print("\n")
        j+=1
        idx=0

#all_vacancies = np.transpose(all_vacancies, (0,2,1,3,4))
k=0

all_start_idxs = all_start_idxs.astype(int)
averages = np.zeros((len(dir_names), iterations, len(times), num_clusters, 3))
#all_vacancies = np.zeros((len(dir_names),iterations, len(times), size, 3))
#print(f"all_vacancies {all_vacancies}")
# averages = np.zeros((len(dir_names), iterations, num_clusters, len(times), 3))


for direc in range(len(all_vacancies)):
    for iteration in range(iterations):
        for idx in range(len(times)):
            for l in range(num_clusters):
                print(f"idx: {all_start_idxs[direc][iteration][idx][(l+1)]}")
                print(np.mean(all_vacancies[direc][iteration][idx][:all_start_idxs[direc][iteration][idx][(l+1)]], axis = 0))
                averages[direc][iteration][idx][l] = np.mean(all_vacancies[direc][iteration][idx][:all_start_idxs[direc][iteration][idx][(l+1)]], axis = 0)

averages = np.transpose(averages, (0,1,4,3,2)) 
print(f"averages.shape: {averages.shape}")
#print(averages[direc][iteration][0][l])

for direc in range(len(averages)):
    all_times_avg = np.mean(all_times[direc], axis=0)
    for iteration in range(iterations):
        for l in range(num_clusters):
            ax.plot(all_times_avg, averages[direc][iteration][0][l], c=lines_colors[k])
            k += 1

ax.plot(all_times_avg, np.ones_like(all_times_avg)*100, linestyle='dashed', c="grey")
ax.set_xlim(0,np.min(np.mean(all_times, axis=1)[:,-1]))
ax.legend(np.arange(1,num_clusters+1).astype(str))
ax.set_xlabel("Times")
ax.set_ylabel("Averge distance (A)")
plt.savefig(f"COM_{dir_names[0]}.png")
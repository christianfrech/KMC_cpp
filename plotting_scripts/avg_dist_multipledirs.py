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

def matano_boltzmann(A, x, t):
    center = 100.5
    lattice_param = 3.5e-10
    D = 5e-11
    vib_freq = 5e-12

    return A/2 * (sci.special.erf((1 - (center - x)) * lattice_param /np.sqrt(2 * math.pi * t * D * vib_freq)) + sci.special.erf((1 + (center - x)) * lattice_param/np.sqrt(2 * math.pi * t * D * vib_freq)))


def parse_line(line):
    toks = line.split(" ")
    lattice_pos = (toks[0])
    x = (toks[1])
    y = (toks[2])
    z = (toks[3])
    return lattice_pos,x,y,z


def read_vac_file(filename, size, dims):
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
            output[idx][0] = (x*2) % (2*dims[0])
            output[idx][1] = (y*2) % (2*dims[1])
            output[idx][2] = (z*2) % (2*dims[2])

        else: 
            output[idx][0] = (x*2 + 1) % (2*dims[0])
            output[idx][1] = (y*2 + 1) % (2*dims[1])
            output[idx][2] = (z*2 + 1) % (2*dims[2])


        idx += 1

    return output

size = 128
iterations = 25
times = np.arange(5,5005,5)
Ts = [300,500,700,900]
dots_colors = ["purple", "blue", "red", "orange", "green"]
lines_colors = ["blue", "red", "orange"]
dims = [16, 16, 12]
idx=0
fig, ax = plt.subplots()
j=0
dir_names = ["fast_and_slow", "pristine_figs"]
data_names = []
all_vacancies = np.zeros((len(dir_names),iterations,len(times),size,3))
all_times = np.zeros((len(dir_names),iterations,len(times)), dtype = float)

rootdir = os.getcwd()
for subdir, dirs, files in os.walk(rootdir):
    for dir in dirs:
        if (("pristine_lattice" in dir) and ("1traj" not in dir)) or ("fast_and_slow" in dir):
            mid_filepath = os.path.join(rootdir, dir)
            mid_filepath = os.path.join(mid_filepath, "vacs")
            print(mid_filepath)
            data_names.append(dir)

            for subdir, dirs, files in os.walk(mid_filepath):
                for file in files:
                    if ("vacancies_output_" in file) and ("rate" not in file):
                        
                        filepath = os.path.join(mid_filepath, file)
                        line_split = file.split("_")
                        iteration = int(line_split[2])
                        step = int(line_split[3])
                        step_time = float(line_split[4])
                        idx = np.where(times == step)[0]
                        #print(f"j: {j}, iteration: {iteration}, idx: {idx}, step: {step}")
                        coords = read_vac_file(filepath, size, dims)
                        all_vacancies[j][iteration][idx] = coords
                        all_times[j][iteration][idx] = step_time
            
            print("\n")
            j+=1
            
all_vacancies = np.transpose(all_vacancies, (0,2,1,3,4))
k=0
print(f"all_vacancies.shape: {all_vacancies.shape}")
for i in range(len(all_vacancies)):
    print(dir_names[i])
    print(f"all_vacancies[i].shape: {all_vacancies[i].shape}")
    all_vacancies_avg = all_vacancies[i].reshape((len(all_vacancies[i]), (len(all_vacancies[i][0])*len(all_vacancies[i][0][0])), len(all_vacancies[i][0][0][0])))
    all_vacancies_avg = np.transpose(all_vacancies_avg, (0,2,1))
    all_times_avg = np.mean(all_times[i], axis=0)
    print(f"all_vacancies_avg: {all_vacancies_avg.shape}")
    averages = []

    for j in range(len(all_vacancies_avg)):
        averages.append(np.mean(all_vacancies_avg[j][2]))

    averages = (np.ones_like(averages)*24 - averages) * (3.502/2) #inverting distance about end of unit cell (z = 15) and scaling by lattice constant (3.502A)
    #print(averages)
    print(f"averages: {averages.shape}")
    print(f"all_times: {all_times.shape}")
    print(f"all_times_avg: {all_times_avg.shape}")
    print(f"i: {i}")
    ax.plot(all_times_avg, averages, c=lines_colors[k])
    print("heatmap gif")
    k += 1

ax.set_xlim(0,np.min(np.mean(all_times, axis=1)[:,-1]))
ax.legend(data_names)
ax.set_xlabel("Times")
ax.set_ylabel("Averge distance (A)")
plt.savefig(f"fastchannel_vs_pristine_avgdist_plot.png")

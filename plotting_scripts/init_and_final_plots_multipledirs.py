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
times = [0, 100, 333.3, 666.67, 1000 ]
Ts = [300,500,700,900]
lines_colors = ["blue", "red", "orange"]
dir_names = ["LiF_figs", "LiO2_figs"]
labels = ["LiF init", "LiF final", "LiO2 init", "LiO2 final"]
shapes = ["+","s"]
times = np.array([5, 5000])
xmax = 0
dims = [16, 16, 12]
idx=0
fig, ax = plt.subplots()
j=0
plot_names = []
dir_names = ["fast_and_slow", "pristine_figs"]
data_names = ["pristine_figs init", "pristine_figs final", "fast_and_slow init", "fast_and_slow final"]
all_vacancies = np.zeros((len(dir_names),iterations,len(times),size,3))
all_times = np.zeros((len(dir_names),iterations,len(times),size,3))


rootdir = os.getcwd()
for subdir, dirs, files in os.walk(rootdir):
    for dir in dirs:
        if (("pristine_lattice" in dir) and ("1traj" not in dir)) or ("fast_and_slow" in dir):
            mid_filepath = os.path.join(rootdir, dir)
            mid_filepath = os.path.join(mid_filepath, "vacs")
            print(mid_filepath)
            
            for subdir, dirs, files in os.walk(mid_filepath):
                for file in files:
                    if ("vacancies_output_" in file) and ("rate" not in file):
                        line_split = file.split("_")
                        step_time = float(line_split[4])
                        step = int(line_split[3])
                        iteration = int(line_split[2])

                        if ((step == times[0]) or (step == times[-1])):
                            line_split = file.split("_")
                            step_time = float(line_split[4])
                            iteration = int(line_split[2])

                            if (step == times[0]):
                                idx = 0
                            else:
                                idx = 1

                            filepath = os.path.join(mid_filepath, file)
                            coords = read_vac_file(filepath, size, dims)
                            all_vacancies[j][iteration][idx] = coords
                            all_times[j][iteration][idx] = step_time
                        
                
            print("\n")
            j+=1


all_vacancies = np.array(all_vacancies)
all_times = np.array(all_times)           
all_vacancies = np.transpose(all_vacancies, (0,2,1,3,4))

for i in range(len(all_vacancies)):
    all_vacancies_avg = all_vacancies[i].reshape((len(all_vacancies[i]), (len(all_vacancies[i][0])*len(all_vacancies[i][0][0])), len(all_vacancies[i][0][0][0])))
    all_vacancies_avg = np.transpose(all_vacancies_avg, (0,2,1))
    
    for j in range(len(all_vacancies_avg)):
        new_bins = np.arange(math.floor(np.min(all_vacancies_avg[j][2])), math.ceil(np.max(all_vacancies_avg[j][2])+2), 1)
        counts, bins = np.histogram(all_vacancies_avg[j][2], bins = new_bins)
        norm_arr = np.ones_like(counts) * size * iterations
        counts_normed = np.divide(counts, norm_arr) * 100
        ax.scatter(new_bins[:-1]*(3.502/2), counts_normed, color = lines_colors[i], marker = shapes[j])
            

xtick_labels = np.arange(0,25,1) * (3.502/2)
xtick_labels = np.round(xtick_labels[::6],0).astype(int)
ax.legend(data_names)
ax.set_xlabel("Averge distance from interface (A)")
ax.set_ylabel("Counts")
ax.set_yscale("log")
ax.set_xticks(xtick_labels)
ax.set_xticklabels(xtick_labels[::-1])
plt.savefig(f"fatchannel_vs_pristine_first_and_final_plot_avg.png")
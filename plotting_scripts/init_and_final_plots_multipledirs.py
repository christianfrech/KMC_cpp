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
    return lattice_pos,x,y,z


def read_vac_file(filename, dims):
    read_file = open(filename, 'r')
    lines = [line for line in read_file]
    size = len(lines)
    output = np.zeros((size, 3))
    idx = 0
    for line in lines:
        lattice_pos,x,y,z = parse_line(line)
        x = int(float(x))
        y = int(float(y))
        z = int(float(z))
        
        if (lattice_pos == "v"):
            output[idx][0] = (x) % (2*dims[0])
            output[idx][1] = (y) % (2*dims[1])
            output[idx][2] = (z) % (2*dims[2])

        else: 
            output[idx][0] = (x) % (2*dims[0])
            output[idx][1] = (y) % (2*dims[1])
            output[idx][2] = (z) % (2*dims[2])


        idx += 1

    return output


def find_max_steps(dir):
    max_time = 0

    for subdir, dirs, files in os.walk(dir):
        for file in files:
            if ("fixed" not in file) and ("00000_" in file): 
                split_filename = file.split("_")
                time = float(split_filename[3])

                if (time > max_time): 
                    max_time = time
    
    return max_time


def find_min_steps(dir):
    min_time = sys.maxsize

    for subdir, dirs, files in os.walk(dir):
        for file in files:
            if ("fixed" not in file): 
                split_filename = file.split("_")
                time = float(split_filename[3])

                if (time < min_time): 
                    min_time = time
    
    return min_time


def find_step_size(dir):
    min_time = sys.maxsize
    prev_min_time = sys.maxsize

    for subdir, dirs, files in os.walk(dir):
        for file in files:
            if ("fixed" not in file): 
                split_filename = file.split("_")
                time = float(split_filename[3])

                if (time < min_time): 
                    prev_min_time = min_time
                    min_time = time
    
    step_size = int(prev_min_time) - int(min_time)

    return step_size



def get_size(filename):
    read_file = open(filename, 'r')
    lines = [line for line in read_file]
    size = len(lines)

    return size

dots_colors = ["purple", "blue", "red", "orange", "green"]
lines_colors = ["blue", "red", "orange", "green"]

dir_names = ["void_bulk_large_3NN_centered"]
labels = ["bulk init", "bulk final", "interface init", "interface final"]
shapes = ["+","s"]

rootdir = os.getcwd()

min_time = find_min_steps((rootdir + "/" + dir_names[0] + "/vacs"))
max_time = find_max_steps((rootdir + "/" + dir_names[0] + "/vacs"))
step_time = find_step_size((rootdir + "/" + dir_names[0] + "/vacs"))
print(f"min_time: {int(min_time)}   max_time: {int(max_time)}   step_time: {int(step_time)}")
times = np.array([int(min_time), int(max_time)])
times_array = np.arange(times[0], times[1]+1, step_time, dtype = int)
dims = [130, 130, 130]
idx = 0
size = 0
fig, ax = plt.subplots()
j = 0
iterations = 1
num_files = int(max_time / step_time) + 1
print(f"num_files: {num_files}")
plot_names = []

for dir in filter(os.path.isdir, os.listdir(os.getcwd())):
    if (dir in dir_names):
        mid_filepath = os.path.join(rootdir, dir)
        mid_filepath = os.path.join(mid_filepath, "vacs")
        print(mid_filepath)

        for subdir, dirs, files in os.walk(mid_filepath):
            for file in files:
                if ("vacancies_output_" in file) and ("rate" not in file) and ("fixed" not in file):
                    filepath = os.path.join(mid_filepath, file)
                    size = get_size(filepath)
                
                break


print(f"size: {size}")

print(f"times: {times}")
all_vacancies = np.zeros((len(dir_names), iterations, len(times), size,  3))
all_times = np.zeros((len(dir_names), iterations, num_files))
print(f"all_vacancies.shape: {all_vacancies.shape}")

for dir in filter(os.path.isdir, os.listdir(os.getcwd())):
    if (dir in dir_names):
        mid_filepath = os.path.join(rootdir, dir)
        mid_filepath = os.path.join(mid_filepath, "vacs")
        print(mid_filepath)

        for subdir, dirs, files in os.walk(mid_filepath):
            for file in files:
                
                if ("vacancies_output_" in file) and ("rate" not in file) and ("fixed" not in file):
                    
                    line_split = file.split("_")
                    step_time = float(line_split[4])
                    step = int(line_split[3])
                    iteration = int(line_split[2])

                    if (step <= times[1]):

                        idx = np.where(times_array == step)[0]
                        if (len(idx) == 0): continue

                        all_times[j][iteration][idx] = step_time
                        idx += 1

        j += 1
        idx = 0

j = 0

for i in range(len(all_times)):  all_times[i][0] = np.sort(all_times[i][0])
maxes = np.array([np.max(all_times[i]) for i in range(len(all_times))])
print(f"maxes: {maxes}")
print(f"np.min(maxes): {np.min(maxes)}")
max_idx = np.where(maxes == np.min(maxes))[0]
print(f"max_idx: {max_idx}")

for i in range(len(all_times)): print(f"np.max(all_times[{i}]): {np.max(all_times[i])}")

print(f"all_times.shape {all_times.shape}")
time_idxs = []

for i in range(len(all_times)):
    if (all_times[i][0][-1] == np.max(all_times[max_idx][0]) ): 
        print(f"appending len: {len(all_times[i][0]) - 1}")
        time_idxs.append(len(all_times[i][0]) - 1) 
    else: time_idxs.append(np.searchsorted(all_times[i][0], np.max(all_times[max_idx][0])))

print(f"time_idxs: {time_idxs}")
print(f"times_array[time_idxs[j]]: {times_array[time_idxs[j]]}")

for dir in filter(os.path.isdir, os.listdir(os.getcwd())):
    if (dir in dir_names):
        plot_names.append(dir.strip("_centered/") + "_init")
        plot_names.append(dir.strip("_centered/") + "_final")
        mid_filepath = os.path.join(rootdir, dir)
        mid_filepath = os.path.join(mid_filepath, "vacs")
        print(mid_filepath)
        print(f"time_idxs[j]: {time_idxs[j]}")
        for subdir, dirs, files in os.walk(mid_filepath):
            for file in files:
                if ("vacancies_output_" in file) and ("rate" not in file) and ("fixed" not in file):
                    line_split = file.split("_")
                    step_time = float(line_split[4])
                    step = int(line_split[3])
                    iteration = int(line_split[2])
                    #print(f"step: {step}")

                    if ((step == times[0])):
                        print("case 1")
                        print(f"file: {file}")
                        print(f"j: {j}, step: {step}")
                        idx = 0
                        filepath = os.path.join(mid_filepath, file)
                        coords = read_vac_file(filepath, dims)
                        print(f"coords.shape: {coords.shape}")
                        all_vacancies[j][iteration][idx] = coords
                        all_times[j][iteration][idx] = step_time
                    
                    elif (step == times_array[time_idxs[j]]):
                        print("case 2")
                        print(f"file: {file}")
                        print(f"j: {j}, step: {step}")
                        idx = 1
                        filepath = os.path.join(mid_filepath, file)
                        coords = read_vac_file(filepath, dims)
                        print(f"coords.shape: {coords.shape}")
                        all_vacancies[j][iteration][idx] = coords
                        all_times[j][iteration][idx] = step_time
        
        print("done\n")
        j+=1

print("transposing")


mid_filepath = os.path.join(rootdir, "void_bulk_large_4NN_centered/")
mid_filepath = os.path.join(mid_filepath, "vacs")
filepath = os.path.join(mid_filepath, "vacancies_output_0_0_0_moves.txt")
coords = read_vac_file(filepath, dims)
all_vacancies[0][0][0] = coords
all_times[0][0][0] = 0

all_vacancies = np.array(all_vacancies)
all_times = np.array(all_times)
all_vacancies = np.transpose(all_vacancies, (0,2,1,3,4))
x_max = 0

for i in range(len(all_vacancies)):
    print("\n")
    print(dir_names[i])
    all_vacancies_avg = all_vacancies[i].reshape((len(all_vacancies[i]), (len(all_vacancies[i][0])*len(all_vacancies[i][0][0])), len(all_vacancies[i][0][0][0])))
    all_vacancies_avg = np.transpose(all_vacancies_avg, (0,2,1))
   
    for j in range(len(all_vacancies_avg)):
        new_bins = np.arange(math.floor(np.min(all_vacancies_avg[j][2])), math.ceil(np.max(all_vacancies_avg[j][2])+2), 2)
        counts, bins = np.histogram(all_vacancies_avg[j][2], bins = new_bins)
        if (new_bins[-1] > x_max): x_max = new_bins[-1]
        norm_arr = np.ones_like(counts) * size * iterations
        counts_normed = np.divide(counts, norm_arr) * 100
        #if (j==0): print(f"all_vacancies_avg[j][2]: {all_vacancies_avg[j][2]}")
        ax.scatter(new_bins[:-1]*(3.502/2), counts_normed, color = lines_colors[j], marker = shapes[j]) #change color = lines_colors[j] back to color = lines_colors[i]

xtick_labels = np.arange(0,x_max+1,1) 
xtick_labels = np.round(xtick_labels[::15]* (3.502/2),0).astype(int)
print(f"xtick_labels: {xtick_labels}")
ax.legend(plot_names, loc='upper left')
ax.set_xlabel("Averge distance from center of void (A)")
ax.set_ylabel("Counts")
ax.set_yscale("log")
ax.text(5,5, f"Time: {maxes[max_idx][0]}s", fontsize = 12)
ax.set_xticks(xtick_labels)
ax.set_xticklabels(xtick_labels[::-1])
#plt.show()
#plt.savefig(f"interface_gb_dft_first_and_final_plot_variedGB.png")
dir_out = dir_names[0].strip("/")
plt.savefig(f"{dir_out}_first_and_final_plot.png", dpi = 95)
import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import math
import scipy as sci
import os
import sys
import glob
import shutil
import fix_cell_edges


def parse_line(line):
    toks = line.split(" ")
    lattice_pos = (toks[0])
    x = (toks[1])
    y = (toks[2])
    z = (toks[3])
    return lattice_pos,x,y,z


def read_vac_file(filename, z_dim):
    read_file = open(filename, 'r')
    lines = [line for line in read_file]
    output = np.zeros((len(lines), 3))

    idx = 0
    for line in lines:
        lattice_pos,x,y,z = parse_line(line)
        lattice_pos = int(lattice_pos)
        x = int(x)
        y = int(y)
        z = int(z)
        
        if (lattice_pos == 0):
            output[idx][0] = (x*2) % (2*z_dim)
            output[idx][1] = (y*2) % (2*z_dim)
            output[idx][2] = (z*2) % (2*z_dim)

        else: 
            output[idx][0] = (x*2 + 1) % (2*z_dim)
            output[idx][1] = (y*2 + 1) % (2*z_dim)
            output[idx][2] = (z*2 + 1) % (2*z_dim)

        #print(f"output: {output[idx]}")
        idx += 1

    return output


def write_expanded_vac_file(filename, lines):
    
    file_out = open(filename, "w")

    file_out.write(f"{len(lines)}\n")
    file_out.write(f"\n")
    
    for line in lines:
        x = line[0]
        y = line[1]
        z = line[2]
        x = str(x)
        y = str(y)
        z = str(z)
        
        file_out.write(f"v {x} {y} {z}\n")

    file_out.close()

a1 = [1,1,-1]
a2 = [-1,1,1]
a3 = [1,-1,1]
pwd = os.getcwd()
rootdir = os.path.join(pwd, "vacs")
dims = [130,130,130]
nums_frames = 1000
frames = np.arange(0,nums_frames,1)
newdir = os.path.join(pwd, "expanded_vacs_bulk_long")
filepaths = []
times = []

if os.path.exists(newdir): 
    if (os.path.exists("expanded_vacs_bulk_long/ovito_vis.py")): shutil.move("expanded_vacs_bulk_long/ovito_vis.py", "./ovito_vis.py")
    files = glob.glob(f"{newdir}/*")
    for f in files: os.remove(f)
    if (os.path.exists("./ovito_vis.py")): shutil.move("./ovito_vis.py", "expanded_vacs_bulk_long/ovito_vis.py")
    
else:
    os.mkdir("expanded_vacs_bulk_long")


for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if ("vacancies_output_" in file) and ("rate" not in file) and ("_0_" in file) and ("00000_" in file):
            line_split = file.split("_")
            times.append(int(line_split[3]))
            filepaths.append(os.path.join(rootdir, file))

idxs = np.argsort(np.array(times))
filepaths = np.array(filepaths)[idxs]

print(f"filepaths: {filepaths}") 
if (len(filepaths) < len(frames)): frames = np.arange(0, len(filepaths), 1)

for i in frames:    
    coords = read_vac_file(filepaths[i], dims[2])
    path = os.path.split(filepaths[i])
    output_extension = path[-1][:-4] + "_expanded.xyz"
    filepath_out = os.path.join(newdir, output_extension)
    print(filepath_out)
    write_expanded_vac_file(filepath_out, coords)

fix_cell_edges.run_fix_edges()

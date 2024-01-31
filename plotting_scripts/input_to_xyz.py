import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import math
import scipy as sci
import os
import sys


def parse_line(line):
    print(line)
    toks = line.split(" ")
    lattice_pos = (toks[0])
    x = (toks[1])
    y = (toks[2])
    z = (toks[3])
    return lattice_pos,x,y,z


def read_vac_file(filename, size, z_dim):
    output = np.zeros((size, 3))
    read_file = open(filename, 'r')
    lines = [line for line in read_file]
    idx = 0
    lines = lines[4:]

    for line in lines:
        lattice_type,x,y,z = parse_line(line)
        
        x = math.floor(float(x))
        y = math.floor(float(y))
        z = math.floor(float(z))
        
        if (lattice_type == "v"):
            output[idx][0] = (x*2)
            output[idx][1] = (y*2)
            output[idx][2] = (z*2)

        else: 
            output[idx][0] = (x*2 + 1)
            output[idx][1] = (y*2 + 1)
            output[idx][2] = (z*2 + 1)

        print(output[idx])
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
        #print(f"{x} {y} {z}")
        file_out.write(f"v {x} {y} {z}\n")

    file_out.close()


pwd = os.getcwd()
rootdir = pwd 
size = 1284
dims = [50,32,64]
newdir = pwd
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if ("1024_vacancies_flat.txt" in file):
            filepath_in = os.path.join(rootdir, file)
            coords = read_vac_file(filepath_in, size, dims[2])
            output_extension = file[:-4] + "_expanded.xyz"
            filepath_out = os.path.join(newdir, output_extension)
            write_expanded_vac_file(filepath_out, coords)

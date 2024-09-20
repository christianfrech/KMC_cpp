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


def read_vac_file(filename):
    read_file = open(filename, 'r')
    lines = [line for line in read_file]
    output = np.zeros((len(lines) - 6, 3))
    idx = 0
    lines = lines[6:]

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
    
    print("done")
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
    print("done writing")
    file_out.close()
    return 0


pwd = os.getcwd()
rootdir = pwd #os.path.join(pwd, "vacs") 
newdir = pwd #os.path.join(pwd, "expanded_vacs")
for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if ("test_hex.txt" in file):
            filepath_in = os.path.join(rootdir, file)
            coords = read_vac_file(filepath_in)
            output_extension = file[:-4] + "_expanded.xyz"
            filepath_out = os.path.join(newdir, output_extension)
            print("writing")
            write_expanded_vac_file(filepath_out, coords)
            break

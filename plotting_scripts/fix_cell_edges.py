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
    toks = line.split(" ")
    lattice_site = toks[0]
    x = float(toks[1])
    y = float(toks[2])
    z = float(toks[3])
    return lattice_site,x,y,z


def read_vac_file(filename, z_dim, corner_sites):
    read_file = open(filename, 'r')
    lines = [line for line in read_file]
    lines = lines[2:]
    output = np.zeros((len(lines), 4))
    idx = 0
    in_corner_sites = np.zeros(8)
    corner_sites_list = np.array(list(corner_sites))
    corner_sites_list = list(tuple(i) for i in corner_sites_list)

    for line in lines:
        
        lattice_pos,x,y,z = parse_line(line)
        x = int(x)
        y = int(y)
        z = int(z)

        if tuple([x,y,z]) in corner_sites: 
            print(filename)
            print(tuple([x,y,z]))
            print(corner_sites_list)
            corner_idx = [tup for tup, val in enumerate(corner_sites_list) if val == tuple([x,y,z])][0]
            print(f"idx: {idx}\n")
            in_corner_sites[corner_idx] = 1
        
        if (lattice_pos == "v"): output[idx][0] = 0
        else: output[idx][0] = 1

        output[idx][1] = x
        output[idx][2] = y
        output[idx][3] = z

        idx += 1
        
    num_pads = len(np.nonzero(np.logical_not(in_corner_sites))[0])
    output = np.pad(output, ((0,num_pads),(0,0)), constant_values=((0, 0),(0, 0)))

    j = 0
    for site in in_corner_sites:
        if (site == 0): 
            output[idx][0] = 1
            output[idx][1] = corner_sites_list[j][0]
            output[idx][2] = corner_sites_list[j][1]
            output[idx][3] = corner_sites_list[j][2]
            idx += 1
        j += 1

    return output


def write_expanded_vac_file(filename, lines):
    filename = filename.strip(".xyz") + "_fixed.xyz"
    file_out = open(filename, "w")
    file_out.write(f"{len(lines)}\n")
    file_out.write(f"\n")
    
    for line in lines:
        if (line[0] == 0): 
            lattice_pos = "v"
            color = "1.0 0.378 0.431 "
        else: 
            lattice_pos = "b"
            color = "1.0 1.0 1.0 "

        x = line[1]
        y = line[2]
        z = line[3]
        x = str(x)
        y = str(y)
        z = str(z)

        file_out.write(f"{lattice_pos} {x} {y} {z} {color}\n")

    file_out.close()


def run_fix_edges():
    pwd = os.getcwd()
    rootdir = os.path.join(pwd, "expanded_vacs_bulk_long")
    size = 4399
    dims = np.array([2*130,2*130,2*130],dtype=int)
    nums_frames = 1000
    frames = np.arange(0,nums_frames,1)
    filepaths = []
    times = []
    corner_sites_list = [[0,0,0],[dims[0],0,0],[dims[0],dims[1],0],[dims[0],dims[1],dims[2]],[0,dims[1],dims[2]],[dims[0],0,dims[2]],[0,dims[1],0],[0,0,dims[2]]]
    corner_sites = set(tuple(i) for i in corner_sites_list)

    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if ("vacancies_output_" in file) and ("rate" not in file) and ("_0_" in file) and ("fixed" not in file): # and ("0000_" in file):
                filepath = os.path.join(rootdir, file)  
                coords = read_vac_file(filepath, dims[2], corner_sites)
                write_expanded_vac_file(filepath, coords)
    
    return 0
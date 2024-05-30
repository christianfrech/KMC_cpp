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
    toks = line.split()
    x = (toks[0])
    y = (toks[1])
    z = (toks[2])
    return x,y,z


def read_poscar(filename):
    print(f"filename: {filename}")
    filepath = os.path.join(os.getcwd(), filename)
    read_file = open(filepath, 'r')
    lines = [line for line in read_file]
    num_atoms = lines[6].split()
    num_atoms = [int(num) for num in num_atoms]
    size = np.sum(num_atoms)
    print(f"size: {size}\n")
    output = np.zeros((size, 3))
    header = lines[:9]
    lines = lines[9:]

    idx = 0
    count = 0

    for line in lines:
        if (count == (size)): break
        x,y,z = parse_line(line)
        x = float(x)
        y = float(y)
        z = float(z)

        output[idx] = [x,y,z]
        idx += 1
        count += 1

    return output, header


def int_to_str(coords):
    coords_str = []
    for coord in coords:
        string_out = ""
        for pos in coord: string_out += f"{(pos):.16f} "

        string_out += " T T T \n"
        coords_str.append(string_out)
    
    print(coords_str)
    return coords_str


for subdir, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if ("POSCAR" in file) and ("centered" not in file):

            read_coords, header = read_poscar(file)
            min_z = np.min(np.transpose(read_coords)[2])
            max_z = np.max(np.transpose(read_coords)[2])
            max_diff = 1 - max_z
            min_diff = min_z
            print(f"max_diff: {max_diff}  min_diff: {min_diff}")
            shift = (max_diff - min_diff) / 2

            shift_arr = shift * np.ones((len(read_coords), 3)) 
            shift_arr[:,:2] = 0
            read_coords = read_coords + shift_arr
            print(f"np.transpose(read_coords)[2]: {np.transpose(read_coords)[2]}")
            print(f"np.max(np.transpose(read_coords)[2]): {np.max(np.transpose(read_coords)[2])}")

            read_coords = int_to_str(read_coords)
            lines = header + read_coords
            filename_components = file.split(".")
            filename_components[0] += "_centered"
            filename = filename_components[0] + "." + filename_components[1]
            print(f"filename: {filename} \n\n")
            file_out = open(filename, "w")
            for line in lines: file_out.write(line)
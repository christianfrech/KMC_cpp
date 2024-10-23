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
    filepath = os.path.join(os.getcwd(), filename)
    read_file = open(filepath, 'r')
    lines = [line for line in read_file]
    num_atoms = lines[6].split()
    num_atoms = [int(num) for num in num_atoms]
    size = np.sum(num_atoms)
    output = np.zeros((size, 3))
    header = lines[:8]
    lines = lines[8:]

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


def add_atoms(coords, header, new_coords, new_ids, old_dims, new_dims):
    atom_ids = header[5].split()
    num_atoms = header[6].split()
    num_atoms = [int(num) for num in num_atoms]
    size = np.sum(num_atoms)
    dict_id = {}

    for i in range(len(atom_ids)):
        id = atom_ids[i]
        num = num_atoms[i]
        dict_id[id] = num

    for i in range(len(coords)):
        coords[i][0] = coords[i][0]/new_dims[0] * old_dims[0]
        coords[i][1] = coords[i][1]/new_dims[1] * old_dims[1]
        coords[i][2] = coords[i][2]/new_dims[2] * old_dims[2]

    for i in range(len(new_coords)):
        coord = new_coords[i]
        id = new_ids[i]
        pos = dict_id[id]
        coords = np.insert(coords, [pos], coord, axis=0)
        print(f"pos: {pos}  coord: {coord}")
        dict_id[id] = dict_id[id] + 1

    
    num_atoms_new = []
    for i in range(len(atom_ids)):
        id = atom_ids[i]
        num_atoms_new.append(dict_id[id])

    num_atoms_string = ""
    for num in num_atoms_new:
        num_atoms_string += str(num)
        num_atoms_string += "  "

    num_atoms_string += "\n"
    header[6] = num_atoms_string

    output_str_list = []
    for coord in coords:
        string_out = ""
        for pos in coord:
            string_out += f"{(pos):.16f}"
            string_out += "  "
        
        string_out += "\n"
        output_str_list.append(string_out)

    return output_str_list, header


def generate_bcc(dims):
    a = 3.502
    a1 = np.array([(a/2),(a/2),-(a/2)])
    a2 = np.array([-(a/2),(a/2),(a/2)])
    a3 = np.array([(a/2),-(a/2),(a/2)])
    dim_a = dims[0]
    dim_b = dims[1]
    dim_c = dims[2]
    coords = np.zeros(((dim_a * dim_b * dim_c), 3))
    idx = 0

    for i in range(dim_a):
        for j in range(dim_b):
            for k in range(dim_c):
                pos = a1 * i + a2 * j + a3 * k
                coords[idx] = pos
                idx += 1

    return coords


shifts = [3.2 / 2, 3.3 / 2, 3.4 / 2, 3.5 / 2, 3.6 / 2, 3.7 / 2, 3.8 / 2]
vec_dims = [10,10,10]
sim_cell_dims = [7.77800, 7.77800, 3]

for subdir, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if ("POSCAR" in file) and ("expanded" not in file):
            for shift in shifts:
                new_coords = generate_bcc(vec_dims)
                coords_transpose = np.transpose(new_coords)
                new_coords = new_coords[(3.502*2 > coords_transpose[0]) & (0 <= coords_transpose[0]) & (3.502*2 > coords_transpose[1]) & (0 <= coords_transpose[1]) & (1 > coords_transpose[2]) & (0 <= coords_transpose[2])]
                print(f"new_coords: {new_coords}")
                new_ids = []
                for i in range(len(new_coords)): new_ids.append("Li")

                top_coords = generate_bcc(vec_dims)
                coords_transpose = np.transpose(top_coords)
                top_coords = top_coords[(3.502*2 > coords_transpose[0]) & (0 <= coords_transpose[0]) & 
                    (3.502*2 > coords_transpose[1]) & (0 <= coords_transpose[1]) & 
                    (sim_cell_dims[2] > coords_transpose[2]) & (1 <= coords_transpose[2])]

                print(f"top_coords: {top_coords}")
                top_ids = []
                for i in range(len(top_coords)): top_ids.append("Li")

                read_coords, header = read_poscar(file)
                z_dim = float((header[4].split())[2]) + 2 * shift

                shift_arr = shift * np.ones((len(read_coords), 3)) / float((header[4].split())[2])
                shift_arr[:,:2] = 0
                read_coords = read_coords + shift_arr
                print(f"np.transpose(read_coords)[2]: {np.transpose(read_coords)[2]}")

                max_z = np.max(np.transpose(read_coords)[2]) * float((header[4].split())[2])
                print(f"np.max(np.transpose(read_coords)[2]): {np.max(np.transpose(read_coords)[2])}")
                print(f"z_dim: {z_dim}")
                print(f"max_z: {max_z}")
                print(f"max_z + shift: {max_z + shift}")
                shift_arr = (max_z) * np.ones((len(new_coords), 3))
                shift_arr[:,:2] = 0
                top_coords = top_coords + shift_arr
                new_coords = np.append(new_coords, top_coords, axis = 0)
                new_ids = new_ids + top_ids

                old_dims = [float((header[2].split())[0]), float((header[3].split())[1]), float((header[4].split())[2])]
                new_dims = [float((header[2].split())[0]), float((header[3].split())[1]), z_dim]

                for i in range(len(new_coords)): 
                    new_coords[i][0] = new_coords[i][0]/new_dims[0]
                    new_coords[i][1] = new_coords[i][1]/new_dims[1]
                    new_coords[i][2] = new_coords[i][2]/new_dims[2]
                
                new_c_line = f"     {0:.16f}    {0:.16f}    {z_dim:.16f}\n"
                header[4] = new_c_line
                new_coords, new_header = add_atoms(read_coords, header, new_coords, new_ids, old_dims, new_dims)
                lines = new_header + new_coords
                filename_components = file.split(".")
                filename_components[0] += "_expanded_"
                filename = filename_components[0] + str(shift) + "." + filename_components[1]
                print(f"filename: {filename} \n\n")
                file_out = open(filename, "w")
                for line in lines: file_out.write(line)

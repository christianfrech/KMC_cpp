import numpy as np 
import matplotlib.pyplot as plt 
import os

dir_num = 0
for dir in filter(os.path.isdir, os.listdir(os.getcwd())):
    if (dir.isnumeric()):
        dir_num += 1

fig, ax = plt.subplots()
all_energies = np.zeros((dir_num), dtype = float)
all_dists = np.zeros((dir_num), dtype = float)
path_dist = 3.36
min_energy = np.inf
max_energy = -np.inf
is_converged = np.zeros((dir_num), dtype = float)

rootdir = os.getcwd()
for dir in filter(os.path.isdir, os.listdir(os.getcwd())):
    if (dir.isnumeric()):
        mid_filepath = os.path.join(rootdir, dir)
        print(mid_filepath)

        for subdir, dirs, files in os.walk(mid_filepath):
            for file in files:
                if ( "OSZICAR" in file ):
                    filepath = os.path.join(mid_filepath, file)
                    file_inst = open(filepath)
                    lines = [line for line in file_inst]
                    lines_rev = reversed(lines)
                    in_loop = True

                    for line in lines_rev:
                        if ("E0" in line): 
                            print(line)
                            splitline = line.split()
                            print(splitline)
                            energy = float(splitline[4])
                            break

                    all_energies[int(dir)] = energy 
                    all_dists[int(dir)] = int(dir) * path_dist
                    file_inst.close()
                
                if ( "stdout" in file):
                    filepath = os.path.join(mid_filepath, file)
                    file_inst = open(filepath)
                    lines = [line for line in file_inst]
                    lines_rev = reversed(lines)
                    in_loop = True

                    for line in lines_rev:
                        if ("reached required accuracy" in line): 
                            is_converged[int(dir)] = 1
                            print(f"{dir} is converged")
                            break


        print("\n")

rootfolder = rootdir.split("/")[-1]
print(f"rootfolder: {rootfolder}")
print(f"all_energies: {all_energies}")
print(f"all_dists: {all_dists}")
print(f"is_converged: {is_converged}")
energies_sorted = [energy for _, energy in sorted(zip(all_energies, all_dists))]
ax.set_xlabel("Distance along path (angstroms)")
ax.set_ylabel("Energy (eV)")
ax.scatter(all_dists, all_energies)
ax.plot(all_dists, all_energies)
plt.tight_layout()
plt.savefig(f"energy_vs_dist_{rootfolder}.png")

import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import math
import scipy as sci
import os
import zipfile
import sys
import shutil


def zipdir(zippath, oldpath):
    zf = zipfile.ZipFile(zippath, "w")
    for dirname, subdirs, files in os.walk(oldpath):
        zf.write(dirname)
        for filename in files:
            zf.write(os.path.join(dirname, filename))
    zf.close()


files_to_copy = ["POSCAR", "KPOINTS", "CONTCAR", "OUTCAR"]
dest = "../binding_energies_100_zip"
num_dirs = [str(i) for i in np.arange(0,12)]
zippeddir = "../binding_energies_100_zip.zip"

if not (os.path.isdir(dest)): os.mkdir(dest) 

for subdir, dirs, files in os.walk(os.getcwd()):
    for dir in dirs:
        if dir in num_dirs:
            
            dirdest = os.path.join(dest, dir)
            if not (os.path.isdir(dirdest)): os.mkdir(dirdest) 

            walkdir = os.path.join(os.getcwd(), dir)
            for subdir, dirs, files in os.walk(walkdir):     
                for file in files:
                    filepath = os.path.join(walkdir, file)
                    filedest = os.path.join(dirdest, file)

                    if (file in files_to_copy):
                        shutil.copyfile(filepath, filedest)
                        print(file)
                        print(filedest)

                    elif (".vasp" in file):
                        shutil.copy(filepath, filedest)
                        print(file)
                        print(filedest)

zipdir(zippeddir, dest)
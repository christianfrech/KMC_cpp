import os 
os.environ['OVITO_GUI_MODE'] = '1' # Request a session with OpenGL support
from ovito.io import import_file 
import math
from ovito.vis import Viewport
import numpy as np

filenames = []
times = []

for subdir, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if ("moves_expanded" in file):
            line_split = file.split("_")
            times.append(int(line_split[3]))
            filenames.append(file)

idxs = np.argsort(np.array(times))
filenames = np.array(filenames)[idxs]
filenames = list(filenames)
print(filenames)

pipeline = import_file(filenames, columns = ['Particle Type', 'Position.X', 'Position.Y', 'Position.Z'])
data_list = []
i = 0
print(pipeline.source.num_frames)

for frame in range(pipeline.source.num_frames):
    print(i)
    data = pipeline.compute(frame)
    data_list.append(data)
    i += 1

vp = Viewport()
vp.type = Viewport.Type.Ortho
vp.camera_pos = (-25, -30, -10)
vp.camera_dir = (2, 3, 1)
vp.fov = math.radians(60.0)
data.cell_[:,0] = (12, 0, 0)  # Cell vector 'a'
data.cell_[:,1] = (0, 12, 0)  # Cell vector 'b'
data.cell_[:,2] = (0, 0, 16)  # Cell vector 'c'
data.cell_[:,3] = (0, 0, 0)     # Cell origin

pipeline.add_to_scene()
vp.render_anim(size=(800,600), filename="LiO2_animation.avi", fps=20)
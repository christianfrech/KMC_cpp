import os 
import numpy as np
# Import everything needed to edit video clips  
from moviepy.editor import *
from decimal import Decimal


def create_textonvid_recurs(base_name, name, times, idxs):
    vid_fps = 20
    print(f"len(idxs): {len(idxs)}")
    print(f"{idxs[0]} {idxs[int(np.floor(len(idxs)/2))]} {idxs[-1]}" )

    if (len(idxs) > 400):
        first_filename = create_textonvid_recurs(base_name, name, times, np.arange(idxs[0], idxs[int(np.floor(len(idxs)/2))], dtype=int))
        second_filename = create_textonvid_recurs(base_name, name, times, np.arange(idxs[int(np.floor(len(idxs)/2))], idxs[-1], dtype=int))

        movie1 = VideoFileClip(first_filename)
        movie2 = VideoFileClip(second_filename)
        movie_out =  concatenate_videoclips([movie1, movie2])
        filename_out = "sim_movie_" + str(times[idxs[0]]) + "_" + str(times[idxs[-1]]) + "._times.mp4"
        movie_out.write_videofile(filename_out, fps=vid_fps)
        movie_out.ipython_display(width = 280)
        os.remove(first_filename)
        os.remove(second_filename)

    else:
        final_clip = 0
        movie = VideoFileClip(base_name) 
        for i in idxs:
            print(i)
            i = int(i)
            # loading video dsa gfg intro video 
            raw_clip = movie.subclip(i/vid_fps,(i+1)/vid_fps)
            # Generate a text clip 
            time = '%.3E' % Decimal(times[i])
            txt_clip = TextClip(f"{time} s", fontsize = 25, color = 'black')  
            # setting position of text in the center and duration will be 10 seconds  
            txt_clip = txt_clip.set_pos("bottom").set_duration(1/vid_fps)
            # Overlay the text clip on the first video clip  
            video = CompositeVideoClip([raw_clip, txt_clip]) 
            # showing video  
            if (final_clip == 0):
                final_clip = video
            else:
                final_clip =  concatenate_videoclips([final_clip, video])

        filename_out = "sim_movie_" + str(times[idxs[0]]) + "_" + str(times[idxs[-1]]) +  "_times.mp4"
        final_clip.write_videofile(filename_out, fps=vid_fps)
        final_clip.ipython_display(width = 280)  

    return filename_out

base_name = "0000"

for subdir, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if (".avi" in file):
            base_name = file
            break
        
final_clip = 0
vid_fps = 20
times_init = []
filenames = []
 
for subdir, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if ("moves_expanded" in file):
            line_split = file.split("_")
            times_init.append(float(line_split[4]))
            filenames.append(file)

times_init.sort()
idxs_init = np.arange(0, len(times_init), dtype=int)
name_init = "sim_movie_" + str(times_init[idxs_init[0]]) + "_" + str(times_init[idxs_init[-1]]) +  "_times.mp4"
create_textonvid_recurs(base_name, name_init, times_init, idxs_init)

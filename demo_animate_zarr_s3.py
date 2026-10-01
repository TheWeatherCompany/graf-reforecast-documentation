"""
demo_animate_zarr_s3.py

s3://twc-graf-reforecast/2026-05-01.canonical/

Loops over every available lead time for one GRAF reforecast case and
variable, plots each one (reusing the reading/gridding/plotting code in
demo_read_zarr_s3.py), and assembles the frames into an
animated GIF at 4 frames per second.  This is meant for eyeballing a
whole forecast's worth of frames at once (e.g. to spot any residual
temporal mis-ordering or other data-quality problems).

Usage:

python demo_animate_zarr_s3.py cyyyymmddhh cvar

where cyyyymmddhh is the YEAR/MONTH/DAY/HOUR of the initial condition,
and cvar is the variable name (see demo_read_zarr_s3.py for
the full list).

Example:
    $python demo_animate_zarr_s3.py 2004010112 t2m

Code by Tom Hamill, The Weather Company, 16 June 2025
"""
import os, sys
import shutil
import tempfile
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

import demo_read_zarr_s3 as gr

FPS = 4
FRAME_DPI = 150 # animation frames don't need the 400-dpi detail of a single still plot

# ------------------------------------------------------------------------

cyyyymmddhh = sys.argv[1] # the desired initial cond. in year-month-day-hour format
cvar = sys.argv[2]

directory_string, max_leadtime = gr.lookup_case_directory(cyyyymmddhh)
if directory_string is None:
    sys.exit()

try:
    gr.valid_variables.index(cvar)
except ValueError:
    gr.print_variable_help(cvar)
    sys.exit()

# ---- plotting info that doesn't vary frame-to-frame (clevs/scale/offset/
#      clabel are the same for every lead time; only the title changes).

plot_info = gr.plotting_information(cyyyymmddhh, '0')
clevs = plot_info['clevs'][cvar]
scale_factor = plot_info['scale_factor'][cvar]
offset_factor = plot_info['offset_factor'][cvar]
clabel = plot_info['clabel'][cvar]

output_gif = 'GRAFR_'+cvar+'_'+cyyyymmddhh+'_animation.gif'

lons, lats = gr.load_static_latlon()

m = gr.make_basemap()
inbox_idx, xg, yg, nn_idx, nn_dist = gr.build_or_load_conus_grid(m, lons, lats)
XG, YG = np.meshgrid(xg, yg)

# ---- Open the zarr store once and reuse it for every frame, rather than
#      reopening/re-decoding Time for each lead time.

zgroup, time_vals_sorted, sort_order = gr.open_case_store(directory_string, cvar)
ntimes = len(time_vals_sorted)
print ('Number of available timesteps (frames) = ', ntimes)

yyyy, mm, dd, hh = gr.splitdate(cyyyymmddhh)
ic_datetime = gr.datetime(yyyy, mm, dd, hh)

# ---- For the 15-min-accumulated precipitation variables, read_timestep()
#      sums 3 consecutive 5-min buckets into a trailing 15-min total.  To
#      actually show that as 15-min-cadence output (rather than a 5-min
#      slideshow where each frame's window overlaps the last by 10 min),
#      step through non-overlapping 3-timestep windows here.

if cvar in gr.precip_accum_vars:
    frame_indices = list(range(3, ntimes, 3))
else:
    frame_indices = list(range(ntimes))

tmpdir = tempfile.mkdtemp(prefix='graf_animate_')
frame_files = []

# ---- Build the figure/axes/colorbar once and reuse them for every frame
#      (see setup_animation_plot()'s docstring) instead of calling
#      plot_data(), which creates a brand-new Figure/Axes each time.

anim_state = gr.setup_animation_plot(m, XG, YG, clevs, cvar, clabel)

try:
    for frame_num, time_index in enumerate(frame_indices):

        valid_time = time_vals_sorted[time_index]
        clead_hours = (valid_time - ic_datetime).total_seconds()/3600.

        data_array = gr.read_timestep(zgroup, sort_order, time_index, cvar)
        gridded = gr.regrid_to_conus(data_array, inbox_idx, nn_idx, nn_dist, xg, yg)

        frame_plot_info = gr.plotting_information(cyyyymmddhh, '%.2f' % clead_hours)
        title = frame_plot_info['title'][cvar]

        frame_file = os.path.join(tmpdir, 'frame_%04d.png' % time_index)
        gr.render_animation_frame(anim_state, (gridded-offset_factor)*scale_factor,
            title, frame_file, dpi=FRAME_DPI)
        frame_files.append(frame_file)

        print ('frame %d/%d, valid_time = %s, lead = %.2f h' % \
            (frame_num+1, len(frame_indices), valid_time, clead_hours))

    print ('Assembling', len(frame_files), 'frames into', output_gif, \
        'at', FPS, 'fps ...')

    # ---- Write frames as independent, full-canvas images rather than via
    #      imageio.mimsave's default GIF encoding, which lets Pillow apply
    #      per-frame diffing/local-palette optimization (each frame patched
    #      against the previous one, using a per-frame-adaptive palette).
    #      That's a known source of viewer-dependent playback glitches
    #      (frame content not resetting/compositing the same way in every
    #      GIF decoder) even when the frames read back correctly here.
    #      optimize=False and disposal=2 force every frame to fully
    #      replace the previous one, independent of any other frame.

    pil_frames = [Image.open(f).convert('RGB') for f in frame_files]
    pil_frames[0].save(
        output_gif,
        format='GIF',
        save_all=True,
        append_images=pil_frames[1:],
        duration=int(round(1000./FPS)),
        loop=0,
        optimize=False,
        disposal=2,
    )
    print ('Done!  Wrote ', output_gif)

finally:
    plt.close(anim_state['fig'])
    shutil.rmtree(tmpdir)

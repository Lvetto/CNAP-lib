from Graphs.FigureBuilder import *
from Graphs.Plots import *
from Graphs.BasicPlots import *
from crystal_maker import *
from Graphs.DataClasses import *
import numpy as np
from old.utils import pddf_calculator
from math import sqrt

# create the figure object
fig = FigureBuilder((1, 3))

# create primitive cells for the 3 cubic lattices
a = 1
cells = [sc_cell(a), bcc_cell(a), fcc_cell(a)]
titles = ["Simple cubic", "Body centered cubic", "Face centered cubic"]

# draw the cubic cells
pointclouds = [PointCloud(np.array(i), size = 1000, alpha=0.8) for i in cells]
plane_data = make_cube_planes(a, 1)
planesets = [SurfaceSet(np.array(plane_data), color=(0.8, 0.8, 0.8), alpha=0.2, edgecolor="k")]
[fig.add_plot(better_particle_plot, ax=fig[i], points=[pointclouds[i]], surfaces=planesets, hide_frames=False, title=title) for i,title in enumerate(titles)]


fig2 =  FigureBuilder((1, 3))

# create a bigger lattice to compute pddf on
lattices = [cubic_lattice_from_cell(a, 3, i) for i in cells]

# compute pddf and compile a list of notable values to mark on the plot
d = [pddf_calculator(0, lattices[i], 0.3*a) for i,_ in enumerate(lattices)]
ps = [[[a, sqrt(2)*a], [100, 100]], [[sqrt(3)*a/2, a], [100, 100]], [[sqrt(2)*a/2, a], [1000, 1000]]]

# draw the bar plot
[fig2.add_plot(bar_plot, ax=fig2[n], bar_data=i[0], crosses=i[2], width=0.3, title=i[1]) for n,i in enumerate(zip(d, titles, ps))]

fig2.show()
 

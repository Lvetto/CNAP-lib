import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.axes import Axes
from Graphs.BasicPlots import *
from matplotlib import cm
from Graphs.DataClasses import *

# A plot object that draws an animation showing a bunch of particles moving around and (part of) their trajectories
class particle_plot(animated_plot, plot_3d):
    def __init__(self, ax, data, size, title= "", traj_frames=0, hide_frames=True):
        #animated_plot.__init__(self, ax, data, title)  # doesn't really do anything currently
        plot_3d.__init__(self, ax, data, title)

        if (hide_frames):
            self.ax.set_axis_off()
        
        self.set_lims()

        # object representing particle positions
        self.scatter = self.ax.scatter(data[0][:,0], data[0][:,1], data[0][:,2], s=size)

        # objects needed to draw trajectories
        self.traj_frames = traj_frames  # lenght of the trajectory in frames
        self.trajectories_data = [[] for _ in range(len(data[0]))]  # keeps track of the last traj_frames positions of a particle
        self.trajectories = [self.ax.plot([], [], [], lw=1)[0] for _ in range(len(data[0]))]    # a list of objects representing the trajecotry lines

    def update(self, frame):
        # Extract data for current frame and convert it into x,y,z arrays
        frame_data = self.data[frame]
        x = frame_data[:, 0]
        y = frame_data[:, 1]
        z = frame_data[:, 2]

        # Plot updated particle positions
        self.scatter._offsets3d = (x, y, z)

        # The rest of this function crashes if we are not saving at least one frame for the trajectories due to the shape of the resulting numpy arrays
        if (self.traj_frames == 0):
            return self.scatter,

        # Keep track of the last traj_frames positions to draw the trajectories
        for i in range(len(self.data[0])):
            self.trajectories_data[i].append(frame_data[i].copy())
            if len(self.trajectories_data[i]) > self.traj_frames:
                self.trajectories_data[i].pop(0)
        
        # Update data in the objects representing traj lines
        for i, trajectory in enumerate(self.trajectories):
            trajectory.set_data(np.array(self.trajectories_data[i])[:, :2].T)
            trajectory.set_3d_properties(np.array(self.trajectories_data[i])[:, 2])
    
        # return updated objects
        return self.scatter, self.trajectories

    def set_lims(self):
        # Set limits as the maximum/minimum values reached in the data for each component
        flattened_data = np.vstack(self.data)

        x_coords = flattened_data[:, 0]
        y_coords = flattened_data[:, 1]
        z_coords = flattened_data[:, 2]

        x_min, x_max = np.min(x_coords), np.max(x_coords)
        y_min, y_max = np.min(y_coords), np.max(y_coords)
        z_min, z_max = np.min(z_coords), np.max(z_coords)

        self.ax.set_xlim(x_min, x_max)
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_zlim(z_min, z_max)
    
    def init_func(self):
        self.trajectories_data = [[] for _ in range(len(self.data[0]))] # clears trajecotry data to avoid issues on reset








class better_particle_plot(plot_3d):
    def __init__(self, ax, points=[], lines=[], surfaces=[], labels=[], title="", hide_frames=True):
        super().__init__(ax, [points, lines, surfaces, labels], title)

        # let the user decide if they want to show the axes
        if (hide_frames):
            self.ax.set_axis_off()
        
        # does nothing at the moment! At some point it should probably be implemented, especially if we want to animate the plot!
        self.set_lims()

        #This version of the code keeps a reference to the drawers for the various elements. It's not strictly necessary, but it becomes so if we want to animate the plots
        self.surfaces = []
        self.draw_surfaces(surfaces)

        self.lines = []
        self.draw_lines(lines)

        self.scatters = []
        self.draw_points(points)

        self.labels = []
        self.draw_labels(labels)

    def draw_points(self, pointclouds):
        for cloud in pointclouds:
            scatter = self.ax.scatter(cloud.xs, cloud.ys, cloud.zs, s=cloud.size, color=cloud.color, alpha=cloud.alpha, marker=cloud.marker, depthshade=cloud.depthshade)
            self.scatters.append(scatter)
    
    def draw_lines(self, linesets):
        for lineset in linesets:
            for xs, ys, zs in zip(lineset.xs, lineset.ys, lineset.zs):
                self.lines.append(self.ax.plot(xs, ys, zs, color=lineset.color, alpha=lineset.alpha, linestyle=lineset.linestyle, linewidth=lineset.linewidth))
        
    def draw_labels(self, labelsets):
        for labelset in labelsets:
            for x,y,z,s in zip(labelset.xs, labelset.ys, labelset.zs, labelset.texts):
                self.labels.append(self.ax.text(x, y, z, s, color=labelset.color, alpha=labelset.alpha, fontname=labelset.fontname, fontsize=labelset.fontsize))
    
    def draw_surfaces(self, surfacesets):
        for surfaceset in surfacesets:
            for x,y,z in zip(surfaceset.xs, surfaceset.ys, surfaceset.zs):
                self.surfaces.append(self.ax.plot_surface(x, y, z, vmin=z.min() * 2, vmax=z.max()*2, color=surfaceset.color, alpha=surfaceset.alpha, edgecolor=surfaceset.edgecolor))

    def init_func(self):
        pass

    def set_lims(self):
        pass


class graph_plot(better_particle_plot):
    def __init__(self, ax, graph, pointsize=100, pointcolor="C0", pointalpha=1, linesize=5, linestyle="-", linecolor="k", linealpha=1, labels=[], title="", hide_frames=True):
        # get data from the graph object
        points = graph.positions
        pairs = graph.unique_bonds
        linedata = np.array([(points[i], points[j]) for i,j in pairs])

        # make points and lines objects
        pointclouds = [PointCloud(points, size=pointsize, color=pointcolor, alpha=pointalpha)]
        linesets = [LineSet(linedata, linewidth=linesize, linestyle=linestyle, color=linecolor, alpha=linealpha)]

        # call the init of better_particle_plot with the data created
        super().__init__(ax, pointclouds, linesets, [], labels, title, hide_frames)




# A plot to show an arbitrary amount of 3d lines
class line_3d_plot(plot_3d):
    def __init__(self, ax, data, title=""):
        super().__init__(ax, data, title)

        self.lines = []

        for line in data:
            x = line[:, 0]
            y = line[:, 1]
            z = line[:, 2]

            self.lines.append(self.ax.plot(x, y, z, lw=1))

class bar_plot(basic_plot):
    def __init__(self, ax, bar_data, crosses=[[], []], width=1, title="", x_title="", y_title=""):
        super().__init__(ax,bar_data, title)
        ax.bar(bar_data[0],bar_data[1], edgecolor="k", width=width)
        ax.plot(crosses[0], crosses[1], "rx")
        ax.set_xlabel(x_title)
        ax.set_ylabel(y_title)



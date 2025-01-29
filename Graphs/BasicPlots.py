import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.axes import Axes

# A basic class implementing the basic functionalities of a plot object
class basic_plot:
    def __init__(self, ax, data, title):
        self.ax = ax
        self.data = data
        self.title = title

        ax.set_title(title)

# A class that implements basic functionalities for animated plots.
# Mostly usefull in the builder class to distinguish between static and animated plots
class animated_plot(basic_plot):
    def __init__(self, ax, data, title):
        super().__init__(ax, data, title)
    
    # this function is called on each frame to update it
    def update(self, frame):
        pass

    # this function is called when the animation starts or resets.
    # It usually clears data to avoid glitches on restart
    def init_func(self):
        pass

# Basic 3d plot funcitonalities. This also handles swapping the axis with a 3d one
class plot_3d(basic_plot):
    def __init__(self, ax, data, title):
        super().__init__(ax, data, title)

        self.ax = self.convert_to_3d(self.ax)
        self.ax.set_title(title)    # needs to be called again because the axis object has been replaced
    
    # 2d and 3d axes are different classes altoghether. To turn one into the other you need to find its info, remove it and substitute it with the correct one with info from the former
    def convert_to_3d(self, ax):
        # extract the info needed from the old axis
        fig = ax.figure
        pos = ax.get_position()

        # sub it with a 3d one with the same info
        ax.remove()
        new_ax = fig.add_axes(pos, projection='3d')

        return new_ax

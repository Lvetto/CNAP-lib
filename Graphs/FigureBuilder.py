import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.axes import Axes
from Graphs.BasicPlots import animated_plot


# Class used to build figures composed of several plots.
# Abstracts away most of the work needed to manage a multiplot figure (supporting animated and regular plots) and provides an easy interface to access/edit each plot
class FigureBuilder:
    def __init__(self, shape, figsize=(10, 6), title=None):

        # These are the objects that do most of the work
        self.fig, self.axs = plt.subplots(*shape, figsize=figsize)

        # This makes sure the interface is consistent for single and multi plot figures
        if isinstance(self.axs, Axes):
            self.axs = np.array([self.axs])

        # Allows users to set a title for the figure
        if title:
            self.fig.suptitle(title)
        
        # A list of the animated plots in the figure. Needed to know which update/reset functions need to be called
        self.anim_list = []

        # A list of all the plots added. Currently only useful to keep a reference to them and save them from the garbage collector
        self.plot_list = []
    
    def anim(self, frames, interval):
        # Runs the animations for every animated plot
        # To show the animation, a reference to the object returned by this function MUST be kept in memory to stop the garbage collector from eating it up!
        return FuncAnimation(self.fig, self.update, frames=frames, interval=interval, init_func=self.init_func)

    def add_plot(self, chart_type, **kwargs):
        # Make sure chart_type is callable to avoid undefined and hard to diagnose behaviour
        if not callable(chart_type):
            raise TypeError(f"chart_type must be a class or at least a callable object, instead is: {type(chart_type)}")
          
        # Try to instance chart_type and handle errors
        try:
            instance = chart_type(**kwargs)
            if isinstance(instance, animated_plot):
                self.anim_list.append(instance)
            else:
                self.plot_list.append(instance)
        except TypeError as e:
            raise TypeError(f"Error while instantiating {chart_type} with args: {kwargs}.\nError: {e}")

    def update(self, frame):
        # Calls the update functions for each animated plot and keeps track of the return values (updated objects)
        out = []
        for i in self.anim_list:
            out += i.update(frame)

        return out
    
    def init_func(self):
        # Calls the reset functions for each animation
        for i in self.anim_list:
            i.init_func()

    def __getitem__(self, item):
        # Allows access to the axs array by using the [] syntax on the object.
        # The if-else is needed to support figures with 1 or multiple rows/columns with the same interface
        if (isinstance(item, tuple)) or (isinstance(item, list)):
            return self.axs[*item]
        else:
            return self.axs[item]
    
    def show(self, frames=0, interval=50):
        if (len(self.anim_list)):
            self.anim_obj = self.anim(frames, interval)

        plt.show()

    @property
    def shape(self):
        # this is a simple property that returns the number of rows and columns in the figure
        return self.axs.shape

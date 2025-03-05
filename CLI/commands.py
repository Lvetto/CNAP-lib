from utils import *
from old.pp_io import read_xyz
from old.utils import adjacency_matrix
from Graphs.FigureBuilder import FigureBuilder
from Graphs.Plots import graph_plot, better_particle_plot
from CNA import *
from matplotlib import pyplot as plt
import os
from random import random
from Graphs.DataClasses import *

class Command:
    """
    Base class for commands, storing metadata and execution logic
    """
    
    def __init__(self, name, description, func, aliases=None):
        """
        Initialize a command
        
        :param name: The primary command name
        :param description: A short description of the command
        :param func: The function to execute
        :param aliases: Alternative names for the command
        """
        self.name = name
        self.description = description
        self.func = func
        self.aliases = aliases or []

    def execute(self, state, args):
        """Execute the command with given arguments"""
        return self.func(state, args)

class Help(Command):
    """
    A class representing the help command
    """
    def __init__(self):
        description = """
                         Return a list of available commands or get info about a specific command.
                         First argument can be empty or an alias for a command  """
        name = "help"
        aliases = ["h"]
        func = self.execute
        super().__init__(name, description, func, aliases)
    
    def execute(self, state, args):
        if args:
            try:
                c = find_from_name(args[0], commands, compile_names(commands))
                print(f"{c.name}, {c.aliases}\n\t{c.description}")
            except:
                print("Couldn't find command")
        else:
            out = f"Available commands:"
            comm_strs = [f"\n{i.name} {[i.name] + i.aliases}\t{i.description}" for i in commands]
            print("\n".join([out] + comm_strs))

class LoadXyz(Command):
    def __init__(self):
        name = "loadXYZ"
        aliases = ["lx", "load"]
        description = """
                         Loads and XYZ file into memory
                         Requires a path to an xyz file as the first argument (from /data)"""
        func = self.load
        super().__init__(name, description, func, aliases)
    
    def load(self, state, args):
        if ("data" not in state.keys()):
            state["data"] = []
        try:
            state["data"].append(read_xyz(fr"data/{args[0]}")[1])
        except IndexError:
            print("Error parsing arguments")
        except FileNotFoundError:
            print("Couldn't locate file")

class addGraphPlot(Command):
    def __init__(self):
        name = "graph_plot"
        aliases = ["gp", "graph"]
        description = """
                        Add a graph plot to a figure
                        First argument is the index for the data to plot (from state[data]) (required)
                        Second is the cutoff distance (required)
                        """
        func = self.make_graph
        super().__init__(name, description, func, aliases)
    
    def make_graph(self, state, args):
        if ("figs" not in state.keys()):
            state["figs"] = []
        state["figs"].append(FigureBuilder((1,1)))
        try:
            data = state["data"][int(args[0])]
            state["figs"][-1].add_plot(graph_plot, (0, 0), graph=Graph(adjacency_matrix(0, data, float(args[1])), node_positions=data), pointsize=300, linesize=1, linealpha=0.6, pointalpha=1)
        except:
            print("Argument error!")

class ShowPlot(Command):
    def __init__(self):
        name = "show_plot"
        aliases = ["show", "sp"]
        description = """
                        Shows the figures in state[figs]
        """
        func = self.show
        super().__init__(name, description, func, aliases)
    
    def show(self, state, args):
        plt.show()

class ReadState(Command):
    def __init__(self):
        name = "read_state"
        aliases = ["rs", "state"]
        description = """
                        Read the value of a state variable
                        First argument is the name of the variable to be read (if empty, a list of all state variables is returned instead)
                        """
        func = self.read_state
        super().__init__(name, description, func, aliases)
    
    def read_state(self, state, args):
        if (not args):
            print(list(state.keys()))
            return
        try:
            name = args[0]
            for n, i in enumerate(state[name]):
                print(f"{n}:\t{i}")
        except:
            print("Argument error")

class ListData(Command):
    def __init__(self):
        name = "list_data"
        aliases = ["ls", "list"]
        description = """
                        List the files available in the data folder
        """
        func = self.list
        super().__init__(name, description, func, aliases)
    
    def list(self, state, args):
        print("Contents of /data/:")
        [print(f"\t\t\t{i}") for i in os.listdir("data/")]

class ComputeFVs(Command):
    def __init__(self):
        name = "compute_fvs"
        aliases = ["fvs"]
        description = """
                        Computes and groups feature vectors for a set of points.
                        First argument is the index of the data set (from state[data]) (required)
                        Second argument is the cutoff distance to use (required)
        """
        func = self.fvs
        super().__init__(name, description, func, aliases)
    
    def fvs(self, state, args):
        try:
            ind = int(args[0])
            cutoff = float(args[1])
        except:
            print("Argument error")
            return
        
        try:
            data = state["data"][ind]
            g = Graph(adjacency_matrix(0, data, cutoff), node_positions=data)
        except:
            print("Error loading the data")
            return
        
        if ("fvs" not in state.keys()):
            state["fvs"] = []
        
        fvs = get_feature_vectors(g)
        keys, grouped = group_fvs(fvs)

        state["fvs"].append((data, keys, grouped))

class MapFVs(Command):
    def __init__(self):
        name = "map_fvs"
        aliases = ["map", "mfvs"]
        description = """
                        Plots a set of points, assigning unique colors to different feature vectors.
                        First argument is an index for the Feature Vectors to plot (from state[fvs]) (required)
        """
        func = self.map_fvs
        super().__init__(name, description, func, aliases)
    
    def map_fvs(self, state, args):
        try:
            data, keys, grouped = state["fvs"][int(args[0])]
        except:
            print("Data error")
            return

        if ("figs" not in state.keys()):
            state["figs"] = []

        state["figs"].append(FigureBuilder((1,1)))

        fig = state["figs"][-1]

        colors = [(random(), random(), random()) for _ in grouped]
        pointclouds = [PointCloud(data[i], size=300, color=j, label=k, alpha=0.8) for i,j,k in zip(grouped, colors, keys)]

        fig.add_plot(better_particle_plot, plot_pos=(0,0), points=pointclouds, hide_frames=False, title="", legend=True)

class ExtractFvs(Command):
    def __init__(self):
        name = "decompose_fvs"
        aliases = ["decompose", "dfvs"]
        description = """
                        Creates a plot showing only Feature Vectors selected by the user
                        First argument is the index of the FVs to study (from state[fvs]) (required)
                        The rest is asked to the user in an interactive dialog
        """
        func = self.extract_fvs
        super().__init__(name, description, func, aliases)
    
    def extract_fvs(self, state, args):
        try:
            data, keys, grouped = state["fvs"][int(args[0])]
        except:
            print("Data error")
            return
        
        if ("figs" not in state.keys()):
            state["figs"] = []

        state["figs"].append(FigureBuilder((1,1)))

        print("Available FVs:")
        [print(f"\t{n}: {i}") for n,i in enumerate(keys)]
        inds = input("\n\nEnter a sequence of indices (separated by whitespaces) to select FVs to add to the plot:\n")

        colors = [(random(), random(), random()) for _ in grouped]
        pointclouds = [PointCloud(data[i], size=300, color=j, label=k, alpha=0.8) for i,j,k in zip(grouped, colors, keys)]
        pointclouds = [i for n, i in enumerate(pointclouds) if str(n) in inds]

        fig = state["figs"][-1]
        fig.add_plot(better_particle_plot, plot_pos=(0,0), points=pointclouds, hide_frames=False, title="", legend=True)

class PlotFv(Command):
    def __init__(self):
        name = "plot_fv"
        aliases = ["pfv"]
        description = """
                        Plots the LAE associated with a Feature Vector
                        First argument is the index of the FV data (from state[fvs]) (required)
                        Second is the cutoff radius (required)
                        The user is then asked which FV(s) to plot
        """
        func = self.plot_fvs
        super().__init__(name, description, func, aliases)

    def plot_fvs(self, state, args):
        try:
            data, keys, grouped = state["fvs"][int(args[0])]
        except:
            print("Data error")
            return
        
        if ("figs" not in state.keys()):
            state["figs"] = []
        
        print("Available FVs:")
        [print(f"\t{n}: {i}") for n,i in enumerate(keys)]
        inds = input("\n\nEnter a sequence of indices (separated by whitespaces) to select FVs to add to the plot:\n")
        inds = [i for i in inds if i != " "]

        state["figs"].append(FigureBuilder((1, len(inds))))
        fig = state["figs"][-1]

        g = Graph(adjacency_matrix(0, data, float(args[1])), node_positions=data)

        for n, i in enumerate(inds):
            i = int(i)
            fig.add_plot(graph_plot, plot_pos=(n, 0), graph=make_neighbors_subgraph(g, grouped[i][0]), pointsize=300, linesize=1, linealpha=0.6, pointalpha=1, title=keys[i])


commands = [Help(), LoadXyz(), addGraphPlot(), ReadState(), ShowPlot(), ListData(), ComputeFVs(), MapFVs(), ExtractFvs(), PlotFv()]


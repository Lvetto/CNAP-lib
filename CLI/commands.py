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
from CNAP_database.Patterns import *
from json import load

class Command:
    """
    Base class for commands, storing metadata and execution logic
    """

    def __init__(self, name, description, aliases=None):
        """
        Initialize a command
        
        :param name: The primary command name
        :param description: A short description of the command
        :param aliases: Alternative names for the command
        """
        self.name = name
        self.description = description
        self.aliases = aliases or []

    def execute(self, state, args):
        pass

    def write_to_state(self, state, key, value):
        if (key not in state.keys()):
            state[key] = [value]
        else:
            state[key].append(value)

    def list_state(self, state, key):
        try:
            if (key=="fvs"):
                for n, i in enumerate(state[key]):
                    print(f"{n}:\t {i[0]} patterns: {i[1]}\n")
            else:
                [print(f"{n}:\t {i}\n") for n,i in enumerate(state[key])]

        except KeyError:
            print(f"No available data in state[{key}]")
            raise ValueError

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
        super().__init__(name, description, aliases)

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
            print("(Most commands have an interactive interface that helps the user choose the correct arguments. To use it, call the command without optional arguments)")
            print("(To exit the program type quit or press ctrl-c)")

class LoadXyz(Command):
    """
    A command to load data from xyz files
    """
    def __init__(self):
        name = "loadXYZ"
        aliases = ["lx", "load"]
        description = """
                         Loads an XYZ file into memory
                         First argument is the name of the file to load (from /data/) (optional)
                         """
        super().__init__(name, description, aliases)

    def execute(self, state, args):
        if (len(args) > 0):
            filepath = args[0]
            filepath = filepath if ".xyz" in filepath else filepath + ".xyz"
        else:
            filenames = [i for i in os.listdir(r"data/") if ".xyz" in i]
            print("Available files:")
            [print(f"{n}:\t {i}") for n,i in enumerate(filenames)]
            ind = int(input("Enter the index of the file you want to read: "))
            filepath = filenames[ind]
        
        try:
            data = read_xyz(fr"data/{filepath}")[1]
            self.write_to_state(state, "data", data)
        except FileNotFoundError:
            print(f"Couldn't open file: /data/{filepath}")

class ListData(Command):
    def __init__(self):
        name = "list_data"
        aliases = ["ls", "list"]
        description = """
                        List the files available in the data folder
        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        print("Contents of /data/:")
        [print(f"\t\t\t{i}") for i in os.listdir("data/")]

class ComputeGraph(Command):
    """ 
    A command to create a graph object from loaded data
    """
    def __init__(self):
        name = "compute_graph"
        aliases = ["cg", "graph"]
        description = """
                        Makes a graph object from coordinates data
                        First argument is the index of the data (from state[data]) (optional)
                        Second argument is the cutoff radius (optional)
                    """
        super().__init__(name, description, aliases)

    def execute(self, state, args):
        if (len(args) > 1):
            ind = int(args[0])
            rc = float(args[1])
        else:
            print("Loaded datasets:")
            self.list_state(state, "data")
            ind = int(input("Enter the index of a dataset: "))
            rc = float(input("Enter a cutoff distance: "))
        try:
            data = state["data"][ind]
            g = Graph(adjacency_matrix(0, data, rc), node_positions=data)
            self.write_to_state(state, "graphs", g)
        except IndexError:
            print("Invalid dataset!")
        except:
            print("Error while processing data!")

class ReadState(Command):
    def __init__(self):
        name = "read_state"
        aliases = ["rs", "state"]
        description = """
                        Read the value of a state variable
                        First argument is the name of the variable to be read (if empty, a list of all state variables is returned instead)
                        """
        super().__init__(name, description, aliases)

    def execute(self, state, args):
        if (not args):
            print(list(state.keys()))
            return
        try:
            name = args[0]
            self.list_state(state, name)
        except:
            print("Argument error")

class MakeGraphPlot(Command):
    """
    A command to plot a graph from loaded data
    """
    def __init__(self):
        name = "graph_plot"
        aliases = ["gp"]
        description = """
                        Add a graph plot to a figure
                        First argument is the index for the graph to plot (from state[graphs]) (optional)
                        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        self.write_to_state(state, "figs", FigureBuilder((1,1)))
        if (args):
            ind = int(args[0])
        else:
            print("Available graphs:")
            self.list_state(state, "graphs")
            ind = int(input("Enter the index of the graph to plot: "))
        try:
            data = state["graphs"][ind]
            state["figs"][-1].add_plot(graph_plot, (0, 0), graph=data, pointsize=300, linesize=1, linealpha=0.6, pointalpha=1)
        except IndexError:
            print("Invalid dataset!")
        except:
            print("Error while processing data")

class ShowPlot(Command):
    def __init__(self):
        name = "show_plot"
        aliases = ["show", "sp"]
        description = """
                        Shows the figures in state[figs]
        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        plt.show()
        plt.close()

class ComputeFVs(Command):
    def __init__(self):
        name = "compute_fvs"
        aliases = ["fvs"]
        description = """
                        Computes and groups feature vectors from a graph
                        First argument is the index of the data set (from state[graphs]) (optional)
        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        if (args):
            ind = int(args[0])
        else:
            print("Available graphs:")
            self.list_state(state, "graphs")
            ind = int(input("Enter an index: "))

        try:
            g = state["graphs"][ind]
            fvs = get_feature_vectors(g)
            keys, grouped = group_fvs(fvs)

            index = read_index("CNAP_database")
            state["CNAP-index"] = index
            indexed_keys = [i["sigs"] for i in state["CNAP-index"]]

            new_keys = keys
            for n, key in enumerate(keys):
                if (str(key) in indexed_keys):
                    t = indexed_keys.index(str(key))
                    name = index[t]["pattern_name"]
                    id = index[t]["pattern_id"]
                    new_keys[n] = f"{id} ({name})"

            keys = new_keys

            
            self.write_to_state(state, "fvs", (g, keys, grouped))
        except IndexError:
            print("Invalid index")
        except:
            print("Error while computing FVs")

class Clear(Command):
    def __init__(self):
        name = "clear"
        aliases = ["cl"]
        description = """
                        Clears the terminal by printing a bunch of whitelines                
        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        [print() for _ in range(100)]

class MapFVs(Command):
    def __init__(self):
        name = "map_fvs"
        aliases = ["map", "mfvs"]
        description = """
                        Plots a set of points, assigning unique colors to different feature vectors.
                        First argument is the index of a set of Feature Vectors (from state[fvs]) (optional)
                        Further arguments are treated as indices for the feature vectors in the set to plot
        """
        super().__init__(name, description, aliases)

    def execute(self, state, args):
        if (len(args) > 1):
            d_ind = int(args[0])
            inds = [int(i) for i in args[1:]]

            try:
                graph, keys, grouped = state["fvs"][d_ind]
            except IndexError:
                print("Invalid index!")
                return

        else:
            if (args):
                d_ind = int(args[0])
            else:
                print("Available FVs sets:")
                self.list_state(state, "fvs")
                d_ind = int(input("Enter the index of the set to plot: "))

            try:
                graph, keys, grouped = state["fvs"][d_ind]
            except IndexError:
                print("Invalid index!")
                return

            print("Available FVs:")
            [print(f"{n}: {i}") for n,i in enumerate(keys)]
            inds = input("Emter the indices of the FVs to plot: ").split(" ")
            if (inds[0] == "a"):
                inds = range(len(keys))
            else:
                inds = [int(i) for i in inds if i != " "]

        inds = np.array(inds)

        self.write_to_state(state, "figs", FigureBuilder((1,1)))
        fig = state["figs"][-1]

        colors = np.array([(random(), random(), random()) for _ in grouped])
        pointclouds = [PointCloud(graph.positions[i], size=300, color=j, label=k, alpha=0.8) for i,j,k in zip(grouped, colors, keys)]

        pointclouds = [i for n,i in enumerate(pointclouds) if n in inds]

        fig.add_plot(better_particle_plot, plot_pos=(0,0), points=pointclouds, hide_frames=False, title="", legend=True)

class PlotFV(Command):
    def __init__(self):
        name = "plot_fv"
        aliases = ["pfv"]
        description = """
                        Plots the LAE associated with a Feature Vector
                        First argument is the index of a set of FVs (from state[fvs]) (optional)
                        The user is then asked which FV(s) to plot
        """
        super().__init__(name, description, aliases)

    def execute(self, state, args):
        if (args):
            ind = int(args[0])
        else:
            self.list_state(state, "fvs")
            ind = int(input("Enter the index of the set to process: "))
        
        try:
            graph, keys, grouped = state["fvs"][ind]
        except:
            print("Data error")
            return

        print("Available FVs:")
        [print(f"\t{n}: {i}") for n,i in enumerate(keys)]
        inds = input("\n\nEnter a sequence of indices (separated by whitespaces) to select FVs to add to the plot: ")
        inds = [int(i) for i in inds if i != " "]
        
        self.write_to_state(state, "figs", FigureBuilder((len(inds), 1)))
        fig = state["figs"][-1]

        for n, i in enumerate(inds):
            fig.add_plot(graph_plot, plot_pos=(n, 0), graph=make_neighbors_subgraph(graph, grouped[i][0]), pointsize=300, linesize=1, linealpha=0.6, pointalpha=1, title=keys[i])

class SplitFVs(Command):
    def __init__(self):
        name = "split_fvs"
        aliases = ["sfv", "split"]
        description = """
                        Plot the local environment and spatial distribution of a set of CNAPs
                        First argument is the index of the set of FVs (from state[fvs]) (optional)
        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        if (args):
            ind = int(args[0])
        else:
            self.list_state(state, "fvs")
            ind = int(input("Enter an index: "))
        
        try:
            graph, keys, grouped = state["fvs"][ind]
        except:
            print("Data error")
            return

        print("Available FVs:")
        [print(f"\t{n}: {i}") for n,i in enumerate(keys)]
        inds = input("\n\nEnter a sequence of indices (separated by whitespaces) to select FVs to add to the plot: ").split(" ")
    
        if (inds[0] == "a"):
            inds = range(len(keys))
        else:
            inds = [int(i) for i in inds if i != " "]

        colors = [(random(), random(), random()) for _ in grouped]

        for i in inds:
            self.write_to_state(state, "figs", FigureBuilder((1, 2), title=keys[i]))
            fig = state["figs"][-1]
            fig.add_plot(graph_plot, plot_pos=0, graph=make_neighbors_subgraph(graph, grouped[i][0]), hide_frames=True, pointsize=100, linesize=1, linealpha=0.6, pointalpha=0.8)
            fig.add_plot(graph_plot, plot_pos=1, graph=graph.make_subgraph(grouped[i]), hide_frames=True, pointsize=50, linesize=1, linealpha=0.6, pointalpha=0.8, pointcolor=colors[i])

class ExportFV(Command):
    def __init__(self):
        name = "export_fv"
        aliases = ["efv", "export"]
        description = """
                        Export a pattern to the CNAP database
                        First argument is the index of the FVs set to choose from (optional)
        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        if (args):
            ind = int(args[0])
        
        else:
            self.list_state(state, "fvs")
            ind = int(input("Enter the index of the set of FVs: "))
        
        try:
            graph, keys, grouped = state["fvs"][ind]
        except:
            print("Error while processing data")
            return
        
        print("Available FVs:")
        [print(f"\t{n}: {i}") for n,i in enumerate(keys)]
        ind = int(input("\n\nEnter an index for the FV to export: "))

        state["CNAP-index"] = read_index("CNAP_database")
        indexed_keys = [i["sigs"] for i in state["CNAP-index"]]

        if (str(keys[ind]) in indexed_keys):
            print("The pattern was already saved!")

        else:
            name = input("Enter a name for the pattern: ")
            #id = input("Enter an id for the pattern: ")
            description = input("Enter a description for the pattern: ")
            sources = input("Enter a list of names for sources this pattern can be found in: ").split(" ")

            ids = [int(i["pattern_id"].replace("CNAP-", "")) for i in state["CNAP-index"]]
            id = max(ids) + 1 if len(ids) > 0 else 0
            id = f"CNAP-{id}"

            t = CNAPPattern(str(keys[ind]), make_neighbors_subgraph(graph, grouped[ind][0]), id, name, description, sources=sources)
            t.save()

            create_index("CNAP_database")

class ShowFromDB(Command):
    def __init__(self):
        name = "show_db"
        aliases = ["sdb"]
        description = """
                        Plot a CNAP from the database
        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        index = read_index("CNAP_database")

        for n, CNAP in enumerate(index):
            print(f"{n}:\t{CNAP["pattern_name"]}:\t{CNAP["pattern_id"]}\n\t\t\t{CNAP["sigs"]}\n")
        
        inds = input("\nEnter a sequence of indices (separated by whitespaces) to select patterns to plot: ")
    
        if (inds[0] == "a"):
            inds = range(len(index))
        else:
            inds = [int(i) for i in inds if i != " "]
        
        for n, i in enumerate(inds):
            self.write_to_state(state, "figs", FigureBuilder((1, 1)))
            fig = state["figs"][-1]

            m_data_path = os.path.join("CNAP_database", index[i]["folder"], "metadata.json")
            with open(m_data_path, "r") as file:
                metadata = load(file)

            id = metadata["pattern_id"]
            name = metadata["pattern_name"]
            desc = metadata["description"]
            sources = metadata["sources"]
            sigs = metadata["sigs"]

            g_data_path = os.path.join("CNAP_database", index[i]["folder"], "graph_data.npz")
            t = np.load(g_data_path)
            g = Graph(t["adj_mat"], node_positions=t["pos"])

            fig.add_plot(graph_plot, plot_pos=(0, 0), graph=g, pointsize=300, linesize=1, linealpha=0.6, pointalpha=1)

            fig.fig.text(0.01, 0.8, f"Id: {id}", fontsize=15)
            fig.fig.text(0.01, 0.75, f"Name: {name}", fontsize=15)
            fig.fig.text(0.01, 0.7, f"Description: {desc}", fontsize=15)
            fig.fig.text(0.01, 0.65, f"Found in: {sources}", fontsize=15)
            fig.fig.text(0.01, 0.65, f"Found in: {sources}", fontsize=15)
            fig.fig.text(0.01, 0.6, f"Signatures: {sigs}", fontsize=15)
            fig.fig.text(0.01, 0.55, f"Number of nodes: {g.number_of_nodes}", fontsize=15)
            fig.fig.text(0.01, 0.5, f"Number of bonds: {g.number_of_unique_bonds}", fontsize=15)

class ShowGraph(Command):
    def __init__(self):
        name = "show_graph"
        aliases = ["sg"]
        description = """
                        Plot a graph
                        First argument is the index of the graph to plot (from state[graphs]) (optional)
        """
        super().__init__(name, description, aliases)
    
    def execute(self, state, args):
        if (args):
            ind = int(args[0])
        else:
            print("Available graphs:")
            self.list_state(state, "graphs")
            ind = int(input("Enter the index of the graph to plot: "))
        
        g = state["graphs"][ind]
        self.write_to_state(state, "figs", FigureBuilder((1,1)))
        fig = state["figs"][-1]

        fig.add_plot(graph_plot, plot_pos=0, graph=g, hide_frames=True, pointsize=50, linesize=1, linealpha=0.6, pointalpha=0.8)

commands = [Help(), LoadXyz(), ListData(), ComputeGraph(), ReadState(), MakeGraphPlot(), ShowPlot(),
            ComputeFVs(), Clear(), MapFVs(), PlotFV(), SplitFVs(), ExportFV(), ShowFromDB(), ShowGraph()]

 
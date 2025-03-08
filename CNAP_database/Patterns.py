import json
import numpy as np
from CNA import Graph
import os

class CNAPPattern:
    def __init__(self, folder_path):
        self.folder = folder_path
        self.load_metadata()
        self.load_graph_data()
    
    def __init__(self, sigs, graph, id, name, description, sources=[], extra="", folder=None):
        self.graph = graph

        self.folder = fr"CNAP_database/{id}" if folder is None else folder

        self.metadata = {}
        self.metadata["sigs"] = sigs
        self.metadata["pattern_id"] = id
        self.metadata["pattern_name"] = name
        self.metadata["description"] = description
        self.metadata["additional_notes"] = extra
        self.metadata["sources"] = sources

    def load_metadata(self):
        with open(fr"{self.folder}/metadata.json", "r") as file:
            self.metadata = json.load(file)

    def load_graph_data(self):
        data = np.load(fr"{self.folder}/graph_data.npz", allow_pickle=True)
        self.adjacency_matrix = data["adj_mat"]
        self.coordinates = data["pos"]
        self.graph = Graph(self.adjacency_matrix, node_positions=self.coordinates)
    
    def save(self):
        # Create the directory if it doesn't exist
        os.makedirs(self.folder, exist_ok=True)

        with open(fr"{self.folder}/metadata.json", "w") as file:
            json.dump(self.metadata, file, indent=4)
        
        self.graph.save(fr"{self.folder}/graph_data.npz")

    def summarize(self):
        print(f"Pattern ID: {self.metadata['pattern_id']}")
        print(f"Name: {self.metadata['pattern_name']}")
        print(f"Description: {self.metadata['description']}")
        #print(f"Datasets observed in: {[ds['dataset_name'] for ds in self.metadata['datasets']]}")     #maybe. Needs some more thought put into it

def create_index(data_root, index_file="index.json", make_gif=False):
    index = []

    for root, dirs, files in os.walk(data_root):
        if "metadata.json" in files:
            metadata_path = os.path.join(root, "metadata.json")

            with open(metadata_path, "r") as file:
                metadata = json.load(file)

            # Extract required fields if present
            entry = {
                "sigs": metadata.get("sigs"),
                "pattern_id": metadata.get("pattern_id"),
                "pattern_name": metadata.get("pattern_name"),
                "folder": os.path.relpath(root, data_root)
            }
            index.append(entry)

            if (make_gif):
                g_data_path = os.path.join(root, "graph_data.npz")
                t = np.load(g_data_path)
                g = Graph(t["adj_mat"], node_positions=t["pos"])

                fps=8
                rot_steps = (0, 4, 0)
                rot0 = (20, 0, 0)

                fig = FigureBuilder((1, 1))
                fig.add_plot(rotating_graph_plot, plot_pos=0, graph=g, rot_steps=rot_steps, rot0=rot0 ,hide_frames=True, pointsize=500, linesize=1, linealpha=0.6, pointalpha=0.8)
                fig.show(frames=90, save=True, path=os.path.join(root, "pattern.gif"), fps=fps)
                del fig


    # Save the index.json file
    index_path = os.path.join(data_root, index_file)
    with open(index_path, "w") as index_f:
        json.dump(index, index_f, indent=2)

def read_index(data_root, index_file="index.json"):
    path = os.path.join(data_root, index_file)

    with open(path, "r") as file:
        t = json.load(file)
    return t

if __name__ == "__main__":
    from Graphs.FigureBuilder import FigureBuilder
    from Graphs.Plots import rotating_graph_plot

    create_index("CNAP_database", make_gif=True)

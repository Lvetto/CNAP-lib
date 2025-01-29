import numpy as np

# Classes used to pass arguments to some plot types. They mostly just need to store the values and do very little computation

class PointCloud:
    def __init__(self, points, size=10, color="C0", alpha=1, edgecolors="face", marker="o", depthshade=True):
        if (len(points) > 0):
            self.xs = points[:,0]
            self.ys = points[:,1]
            self.zs = points[:,2]
        else:
            self.xs = []
            self.ys = []
            self.zs = []

        self.size = size
        self.color = color
        self.alpha = alpha
        self.edgecolors = edgecolors
        self.marker = marker
        self.depthshade = depthshade

class LineSet:
    def __init__(self, lines, color="k", alpha=1, linestyle="-", linewidth=1):
        self.xs = [line[:,0] for line in lines]
        self.ys = [line[:,1] for line in lines]
        self.zs = [line[:,2] for line in lines]

        self.color = color
        self.alpha = alpha
        self.linestyle = linestyle
        self.linewidth = linewidth

class LabelSet:
    def __init__(self, positions, texts, fontname=None, fontsize="medium", color="k", alpha=1):
        self.xs = positions[:,0]
        self.ys = positions[:,1]
        self.zs = positions[:,2]

        self.texts = texts
        self.fontname = fontname
        self.fontsize = fontsize
        self.color = color
        self.alpha = alpha

class SurfaceSet:
    def __init__(self, surfaces, color="C0", alpha=None, edgecolor=None):
        self.xs = surfaces[:,0]
        self.ys = surfaces[:,1]
        self.zs = surfaces[:,2]

        self.color = color
        self.alpha = alpha
        self.edgecolor = edgecolor


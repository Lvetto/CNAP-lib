import numpy as np

def write_xyz(particles, filepath):
    with open(filepath, "w") as file:
        file.write(f"\t{len(particles)}\t\n")
        file.write(f"\tFile generato da crystal maker\t\n")
        file.writelines([f"\tX\t{i[0]:.6f}\t{i[1]:.6f}\t{i[2]:.6f}\n" for i in particles])

def cubic_lattice_from_cell(step, repetitions, cell):
    vertices = []

    for x in range(repetitions):
        for y in range(repetitions):
            for z in range(repetitions):
                for v in cell:
                    vertices.append(np.array((x,y,z))*step+v)

    unique_verts = list({tuple(array) for array in (arr.tolist() for arr in vertices)})

    return np.array(unique_verts)

def sc_cell(side_lenght):
    vertices = []

    # basic cube vertices
    vertices.append(np.array((0, 0, 0)))
    vertices.append(np.array((side_lenght, 0, 0)))
    vertices.append(np.array((side_lenght, side_lenght, 0)))
    vertices.append(np.array((0, side_lenght, 0)))
    vertices.append(np.array((0, 0, side_lenght)))
    vertices.append(np.array((side_lenght, 0, side_lenght)))
    vertices.append(np.array((side_lenght, side_lenght, side_lenght)))
    vertices.append(np.array((0, side_lenght, side_lenght)))

    return vertices

def bcc_cell(side_lenght):
    vertices = []

    # basic cube vertices
    vertices.append(np.array((0, 0, 0)))
    vertices.append(np.array((side_lenght, 0, 0)))
    vertices.append(np.array((side_lenght, side_lenght, 0)))
    vertices.append(np.array((0, side_lenght, 0)))
    vertices.append(np.array((0, 0, side_lenght)))
    vertices.append(np.array((side_lenght, 0, side_lenght)))
    vertices.append(np.array((side_lenght, side_lenght, side_lenght)))
    vertices.append(np.array((0, side_lenght, side_lenght)))
    vertices.append(np.array((side_lenght/2, side_lenght/2, side_lenght/2)))

    return vertices

def fcc_cell(side_lenght):
    vertices = []

    # basic cube vertices
    vertices.append(np.array((0, 0, 0)))
    vertices.append(np.array((side_lenght, 0, 0)))
    vertices.append(np.array((side_lenght, side_lenght, 0)))
    vertices.append(np.array((0, side_lenght, 0)))
    vertices.append(np.array((0, 0, side_lenght)))
    vertices.append(np.array((side_lenght, 0, side_lenght)))
    vertices.append(np.array((side_lenght, side_lenght, side_lenght)))
    vertices.append(np.array((0, side_lenght, side_lenght)))

    vertices.append(np.array((side_lenght/2, 0, side_lenght/2)))
    vertices.append(np.array((side_lenght/2, side_lenght, side_lenght/2)))
    vertices.append(np.array((side_lenght/2, side_lenght/2, 0)))
    vertices.append(np.array((side_lenght/2, side_lenght/2, side_lenght)))
    vertices.append(np.array((0, side_lenght/2, side_lenght/2)))
    vertices.append(np.array((side_lenght, side_lenght/2, side_lenght/2)))

    return vertices

def make_cube_planes(step, repetitions):
    """# testing plane drawing functionality
    c = 5
    x = np.linspace(-10, 10, 2)
    y = np.linspace(-10, 10, 2)
    X, Y = np.meshgrid(x, y)
    Z = np.full_like(X, c)"""
    planes = []

    # x-y planes
    for n in range(repetitions+1):
        z = step * n
        x = np.linspace(0, step*repetitions, 2)
        y = np.linspace(0, step*repetitions, 2)
        X, Y = np.meshgrid(x, y)
        Z = np.full_like(X, z)

        planes.append((X,Y,Z))

    # x-z planes
    for n in range(repetitions+1):
        y = step * n
        x = np.linspace(0, step*repetitions, 2)
        z = np.linspace(0, step*repetitions, 2)
        X, Z = np.meshgrid(x, z)
        Y = np.full_like(X, y)

        planes.append((X,Y,Z))
    
    # y-z planes
    for n in range(repetitions+1):
        x = step * n
        y = np.linspace(0, step*repetitions, 2)
        z = np.linspace(0, step*repetitions, 2)
        Y, Z = np.meshgrid(y, z)
        X = np.full_like(Y, x)

        planes.append((X,Y,Z))
    
    return planes


if __name__ == "__main__":
    basic_cell = fcc_cell(1)
    particles = cubic_lattice_from_cell(1, 1, basic_cell)
    write_xyz(particles, "out.xyz")


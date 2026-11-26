"""A space truss tower (a tetrahelix) and one force polyhedron per node.

Form diagram: a tower of tetrahedra on three pinned supports. Every new node is braced to the
three nodes below it, so the truss is statically determinate (3 equations per free node, 3 new
members per free node). Every free node carries gravity plus a horizontal wind load.

Force diagram: the "disjointed" 3GS approach. Each node gets its own closed polyhedral cell, with
one face per force at that node (members, load or reaction). A member appears in the cells of
both its nodes; equilibrium makes those two faces equal in area and parallel. Putting all node
cells together is what the 3GS cell network does.

Steps:
1. Solve the member forces with the equilibrium matrix (numpy), and the support reactions.
2. For every node with at least four forces that span 3D: EGI -> cell -> planarise to areas.
3. Check every face (area = force, normal parallel to force) and every member (same area in the
   cells of both its nodes).

Run with the compas-3gs env:  python examples/physical_structures/20_space_truss.py

"""

import numpy as np
from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas.colors import Color
from compas.geometry import Frame
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Scale
from compas.geometry import Translation
from compas.geometry import add_vectors
from compas.geometry import angle_vectors
from compas.geometry import centroid_points
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas_3gs.algorithms import cell_from_egi
from compas_3gs.algorithms import cell_planarise
from compas_3gs.algorithms import egi_from_vectors

# ------------------------------------------------------------------------------
# 1. geometry: three supports and a tetrahelix of free nodes
# ------------------------------------------------------------------------------

LEVELS = 6  # free nodes
GRAVITY = -2.0  # kN per free node
WIND = 0.5  # kN per free node, along x
KMAX = 2000  # nodes with six or more forces stall at about 1-4 % area error whatever kmax is
LOAD_KEY = 1000  # egi_from_vectors needs integer keys (it numbers new vertices from the largest key)
REACTION_KEY = 2000

supports = [0, 1, 2]
xyz = [[2.0 * np.cos(a), 2.0 * np.sin(a), 0.0] for a in np.radians([90, 210, 330])]
side = np.linalg.norm(np.subtract(xyz[0], xyz[1]))
xyz.append(list(np.mean(xyz, axis=0) + [0, 0, side * np.sqrt(2.0 / 3.0)]))  # regular tetrahedron
for k in range(4, 3 + LEVELS):
    # the next node is the reflection of node k-4 through the face (k-3, k-2, k-1)
    xyz.append(list(2.0 / 3.0 * (np.array(xyz[k - 1]) + xyz[k - 2] + xyz[k - 3]) - xyz[k - 4]))

free = list(range(3, 3 + LEVELS))
members = [(k, k - i) for k in free for i in (1, 2, 3)]
loads = {k: [WIND, 0.0, GRAVITY] for k in free}

# ------------------------------------------------------------------------------
# 2. statics: equilibrium matrix, tension positive
# ------------------------------------------------------------------------------

row = {node: 3 * i for i, node in enumerate(free)}
A = np.zeros((3 * len(free), len(members)))
for m, (i, j) in enumerate(members):
    u = normalize_vector(subtract_vectors(xyz[j], xyz[i]))  # direction from i to j
    if i in row:
        A[row[i] : row[i] + 3, m] += u
    if j in row:
        A[row[j] : row[j] + 3, m] -= u
b = -np.concatenate([loads[k] for k in free])
N = np.linalg.solve(A, b)
print("statically determinate: {} equations, {} members, residual {:.1e}".format(A.shape[0], A.shape[1], np.linalg.norm(A @ N - b)))

# the forces acting on every node, keyed by member index, plus the load or the reaction
forces = {node: {} for node in range(len(xyz))}
for m, (i, j) in enumerate(members):
    u = normalize_vector(subtract_vectors(xyz[j], xyz[i]))
    forces[i][m] = scale_vector(u, float(N[m]))
    forces[j][m] = scale_vector(u, -float(N[m]))
for k in free:
    forces[k][LOAD_KEY] = loads[k]
for s in supports:
    forces[s][REACTION_KEY] = scale_vector([sum(f[c] for f in forces[s].values()) for c in range(3)], -1.0)

print("\nmember forces (kN, + tension / - compression):")
for m, (i, j) in enumerate(members):
    print("  {:2d}: {}-{}  {:+.3f}".format(m, i, j, N[m]))
print("\nsupport reactions (kN):")
for s in supports:
    print("  node {}: {}".format(s, [round(float(c), 3) for c in forces[s][REACTION_KEY]]))
total_load = np.sum([loads[k] for k in free], axis=0)
total_reaction = np.sum([forces[s][REACTION_KEY] for s in supports], axis=0)
print("global equilibrium: loads + reactions = {}".format(np.round(total_load + total_reaction, 9)))

# ------------------------------------------------------------------------------
# 3. one force polyhedron per node
# ------------------------------------------------------------------------------

cells = {}
for node, node_forces in forces.items():
    vectors = {key: f for key, f in node_forces.items() if length_vector(f) > 1e-9}
    try:
        egi = egi_from_vectors(vectors, xyz[node])
    except Exception as error:
        print("node {}: no force polyhedron ({} forces: {})".format(node, len(vectors), error))
        continue
    cell = cell_from_egi(egi)
    targets = {key: length_vector(f) for key, f in vectors.items()}
    for face in cell.faces():
        targets.setdefault(face, 0.0)
    cell_planarise(cell, kmax=KMAX, target_areas=targets, target_normals={key: egi.vertex_attribute(key, "normal") for key in cell.faces()}, tolerance_area=1e-4)
    cells[node] = (cell, vectors)

print("\nforce polyhedra:")
worst_area = worst_angle = 0.0
for node, (cell, vectors) in cells.items():
    errors = [abs(cell.face_area(key) - length_vector(f)) / length_vector(f) for key, f in vectors.items()]
    angles = [angle_vectors(cell.face_normal(key), f, deg=True) for key, f in vectors.items()]
    worst_area, worst_angle = max(worst_area, max(errors)), max(worst_angle, max(angles))
    print("  node {}: {} forces, {} faces, max area error {:.2%}, max angle {:.3f} deg".format(node, len(vectors), cell.number_of_faces(), max(errors), max(angles)))

print("\nmembers seen from both ends (face areas in the two node cells):")
for m, (i, j) in enumerate(members):
    if i in cells and j in cells:
        print("  {:2d}: |N| {:.3f}   cell {}: {:.3f}   cell {}: {:.3f}".format(m, abs(N[m]), i, cells[i][0].face_area(m), j, cells[j][0].face_area(m)))

# ------------------------------------------------------------------------------
# 4. viewer: the tower (blue compression, red tension) and the node cells beside it
# ------------------------------------------------------------------------------

viewer = Viewer()
peak = max(abs(N))
for m, (i, j) in enumerate(members):
    color = Color.blue() if N[m] < 0 else Color.red()
    viewer.scene.add(Line(xyz[i], xyz[j]), linecolor=color, linewidth=1 + 5 * abs(N[m]) / peak)
    viewer.scene.add(Tag("{:+.2f}".format(N[m]), centroid_points([xyz[i], xyz[j]]), color=color, height=12))
for s in supports:
    viewer.scene.add(Point(*xyz[s]), pointcolor=Color.grey(), pointsize=14)
for k in free:
    viewer.scene.add(Line(xyz[k], add_vectors(xyz[k], loads[k])), linecolor=Color.black(), linewidth=2)

SHIFT = [8.0, 0.0, 0.0]
SIZE = 0.35  # cells drawn at 35 % so they do not overlap
for node, (cell, vectors) in cells.items():
    centre = centroid_points([cell.vertex_coordinates(v) for v in cell.vertices()])
    cell.transform(Translation.from_vector(add_vectors(subtract_vectors(xyz[node], centre), SHIFT)) * Scale.from_factors([SIZE] * 3, frame=Frame(centre, [1, 0, 0], [0, 1, 0])))
    facecolor = {}
    for key in cell.faces():
        if key not in vectors:
            facecolor[key] = Color.white()  # faces from arc crossings, no force
        elif key in (LOAD_KEY, REACTION_KEY):
            facecolor[key] = Color.black()
        else:
            facecolor[key] = Color.blue() if N[key] < 0 else Color.red()
    viewer.scene.add(cell, facecolor=facecolor, opacity=0.6, show_lines=True)

viewer.show()

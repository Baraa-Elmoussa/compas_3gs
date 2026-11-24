"""A tripod: three legs carrying a load at the apex, and its force polyhedron.

Form diagram: the apex node, three legs to pinned feet, and an applied load.
Force diagram: one closed polyhedral cell for the apex node. Each face belongs to one force at
the node: its normal is the direction of the force, its area is the magnitude. Because the forces
are in equilibrium, the cell closes (the area-weighted face normals sum to zero).

Steps:
1. Solve the apex equilibrium for the three leg forces (3 equations, 3 unknowns).
2. Give the four force vectors (3 legs + load) to 3GS: EGI -> unit cell -> planarise to areas.
3. Check: every face area equals its force, every face normal is parallel to its force (so each
   face is perpendicular to its member: the reciprocal relation of graphic statics).

Run with the compas-3gs env:  python examples/physical_structures/10_tripod.py

"""

import numpy as np
from compas.colors import Color
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Translation
from compas.geometry import add_vectors
from compas.geometry import angle_vectors
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas_3gs.algorithms import cell_from_egi
from compas_3gs.algorithms import cell_planarise
from compas_3gs.algorithms import egi_from_vectors

# ------------------------------------------------------------------------------
# 1. structure and statics (tension positive)
# ------------------------------------------------------------------------------

HEIGHT = 3.0
RADIUS = 2.0
LOAD = [1.5, 0.0, -10.0]  # kN: gravity plus a horizontal push along x

apex = [0.0, 0.0, HEIGHT]
feet = [[RADIUS * np.cos(a), RADIUS * np.sin(a), 0.0] for a in np.radians([0, 120, 240])]
directions = [normalize_vector(subtract_vectors(foot, apex)) for foot in feet]

# sum_i N_i * u_i + P = 0
N = np.linalg.solve(np.array(directions).T, -np.array(LOAD))

print("leg forces (kN, + tension / - compression):")
for i, n in enumerate(N):
    print("  leg {}: {:+.3f}".format(i, n))

# the forces acting on the apex node: the members and the load
forces = {i: scale_vector(directions[i], float(N[i])) for i in range(3)}
forces[3] = LOAD
print("sum of forces at the apex: {:.2e}".format(length_vector([sum(f[k] for f in forces.values()) for k in range(3)])))

# ------------------------------------------------------------------------------
# 2. force polyhedron with 3GS
# ------------------------------------------------------------------------------

egi = egi_from_vectors(forces, apex)
cell = cell_from_egi(egi)
target_areas = {key: length_vector(force) for key, force in forces.items()}
for face in cell.faces():
    target_areas.setdefault(face, 0.0)  # faces from arc crossings, if any, have no force
cell_planarise(cell, kmax=3000, target_areas=target_areas, target_normals={key: egi.vertex_attribute(key, "normal") for key in cell.faces()})

# ------------------------------------------------------------------------------
# 3. check: area = magnitude, normal parallel to force
# ------------------------------------------------------------------------------

print("\nforce polyhedron of the apex ({} faces):".format(cell.number_of_faces()))
names = {0: "leg 0", 1: "leg 1", 2: "leg 2", 3: "load"}
for key, force in forces.items():
    area = cell.face_area(key)
    angle = angle_vectors(cell.face_normal(key), force, deg=True)
    print("  {:6s} |F| = {:7.3f}   face area = {:7.3f}   angle(normal, F) = {:.2f} deg".format(names[key], length_vector(force), area, angle))

# ------------------------------------------------------------------------------
# 4. viewer: form diagram on the left, force polyhedron on the right
# ------------------------------------------------------------------------------

viewer = Viewer()
shift = [6.0, 0.0, 0.0]
color = {0: Color.red(), 1: Color.green(), 2: Color.blue(), 3: Color.black()}

for i, foot in enumerate(feet):
    viewer.scene.add(Line(apex, foot), linecolor=color[i], linewidth=4)
    viewer.scene.add(Point(*foot), pointcolor=Color.grey(), pointsize=12)
    viewer.scene.add(Tag("{:+.2f} kN".format(N[i]), [(a + b) / 2 for a, b in zip(apex, foot)], color=color[i], height=18))
viewer.scene.add(Line(add_vectors(apex, scale_vector(LOAD, -0.2)), apex), linecolor=Color.black(), linewidth=3)
viewer.scene.add(Tag("P", add_vectors(apex, scale_vector(LOAD, -0.2)), color=Color.black(), height=18))

cell.transform(Translation.from_vector(subtract_vectors(shift, apex)))
viewer.scene.add(cell, facecolor={key: color.get(key, Color.white()) for key in cell.faces()}, opacity=0.7, show_lines=True)
for key in forces:
    viewer.scene.add(Tag("{:.2f}".format(cell.face_area(key)), cell.face_centroid(key), color=Color.black(), height=16))

viewer.show()

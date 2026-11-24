"""Modify a polyhedral force diagram: set vertex fixities, move a vertex and lift a boundary vertex.

A move keeps the fixed axes of the vertex: the moved bottom vertex stays in z = 0, so the bottom
stays flat (the vertical faces through it are warped). The lift adds a pyramid cell (red) on each
boundary halfface around the vertex. Every fixed vertex is labelled with its fixed axes, in red
for x, blue for z and magenta for both. The grey point is the requested position of the moved
vertex, and the arrow its constrained move.
"""

import json
import os

from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas.colors import Color
from compas.geometry import Point
from compas.geometry import Vector
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.geometry import subtract_vectors
from compas_3gs import DATA
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.operations import volmesh_vertex_lift

X_FIXED_VKEYS = [0, 3, 6, 7, 10, 11, 24, 25, 26]  # the vertices in the side x = 0
Z_FIXED_VKEYS = [0, 1, 4, 6, 8, 10, 12, 14, 16]  # the vertices in the bottom z = 0
MOVE_VKEY = 4  # the centre of the bottom, z fixed
MOVE_XYZ = [7.0, 7.0, -3.0]
LIFT_VKEY = 21  # the centre of the top, on the boundary
LIFT_XYZ = [5.0, 5.0, 15.0]  # on the vertex normal

# ------------------------------------------------------------------------------
#   1. force diagram from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_cubes_4x4.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

# ------------------------------------------------------------------------------
#   2. fixity
# ------------------------------------------------------------------------------
for vkey in X_FIXED_VKEYS:
    forcediagram.vertex_attribute(vkey, "x_fix", True)

for vkey in Z_FIXED_VKEYS:
    forcediagram.vertex_attribute(vkey, "z_fix", True)

# ------------------------------------------------------------------------------
#   3. move a vertex (constrained: its fixed axes do not change)
# ------------------------------------------------------------------------------
before = forcediagram.vertex_coordinates(MOVE_VKEY)
forcediagram.vertex_update_xyz(MOVE_VKEY, MOVE_XYZ)
after = forcediagram.vertex_coordinates(MOVE_VKEY)

print("move vertex", MOVE_VKEY, "to", MOVE_XYZ, ":", before, "->", after)

# ------------------------------------------------------------------------------
#   4. lift a boundary vertex, with a new cell on each of its boundary halffaces
# ------------------------------------------------------------------------------
hfkeys = [hfkey for hfkey in forcediagram.vertex_halffaces(LIFT_VKEY) if forcediagram.is_halfface_on_boundary(hfkey)]
old_cells = list(forcediagram.cells())

volmesh_vertex_lift(forcediagram, LIFT_VKEY, LIFT_XYZ, hfkeys)

new_cells = [cell for cell in forcediagram.cells() if cell not in old_cells]

print("lift vertex", LIFT_VKEY, "to", LIFT_XYZ, "on halffaces", hfkeys, ":", len(old_cells), "->", forcediagram.number_of_cells(), "cells")

# ------------------------------------------------------------------------------
#   5. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# cells, the new ones in red
for cell in forcediagram.cells():
    color = Color.red().lightened(50) if cell in new_cells else Color.white()
    viewer.scene.add(forcediagram.cell_to_mesh(cell), facecolor=color, show_lines=True, opacity=0.5, name="cell {}".format(cell))

# vertices, with their fixed axes as labels (x red, y green, z blue, mixed when several are fixed)
for vkey in forcediagram.vertices():
    xyz = forcediagram.vertex_coordinates(vkey)
    viewer.scene.add(Point(*xyz), pointcolor=Color.black(), pointsize=10)
    text = ""
    rgb = [0, 0, 0]
    for i, axis in enumerate("xyz"):
        if forcediagram.vertex_attribute(vkey, axis + "_fix"):
            text += axis
            rgb[i] = 255
    if text:
        viewer.scene.add(Tag(text, xyz, color=Color.from_rgb255(*rgb), height=20))

# the move: the requested target (grey) and the constrained move (arrow)
viewer.scene.add(Point(*MOVE_XYZ), pointcolor=Color.grey(), pointsize=10)
viewer.scene.add(Vector(*subtract_vectors(after, before)), anchor=Point(*before), linecolor=Color.black())

# camera on the diagram and the requested position (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz") + [MOVE_XYZ])
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10

viewer.show()

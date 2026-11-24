"""Construct a force volmesh from a set of cells and show its cells, vertices and vertex keys."""

import json
import os

from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas.colors import Color
from compas.geometry import Point
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas_3gs import DATA
from compas_3gs.diagrams import ForceVolMesh

# ------------------------------------------------------------------------------
#   1. force volmesh from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_cubes_4x4.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

print("cells:", forcediagram.number_of_cells())
print("vertices:", forcediagram.number_of_vertices())
print("halffaces:", len(list(forcediagram.halffaces())))

# ------------------------------------------------------------------------------
#   2. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# cells: one mesh per cell, white faces
for cell in forcediagram.cells():
    mesh = forcediagram.cell_to_mesh(cell)
    viewer.scene.add(mesh, facecolor={fkey: Color.white() for fkey in mesh.faces()}, show_lines=True, opacity=0.5, name="cell {}".format(cell))

# vertices: points with their keys as tags
for vkey in forcediagram.vertices():
    xyz = forcediagram.vertex_coordinates(vkey)
    viewer.scene.add(Point(*xyz), pointcolor=Color.black(), pointsize=10)
    viewer.scene.add(Tag(str(vkey), xyz, color=Color.black(), height=20))

# camera on the diagram (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10

viewer.show()

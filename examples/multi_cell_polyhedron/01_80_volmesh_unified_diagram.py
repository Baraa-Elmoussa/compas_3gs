"""Unified diagram of a multi-cell force diagram and its form diagram.

The form diagram is reciprocated to the force diagram, then every cell is scaled towards its
form node and neighbouring cells are joined by prisms. The cells are black; the prisms and the
form edges are coloured by force magnitude (blue compression, red tension, darker is larger).
The smallest force is drawn white, so the background is grey.
"""

import json
import os

from compas_viewer import Viewer

from compas.colors import Color
from compas.datastructures import Mesh
from compas.geometry import Line
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.algorithms import volmesh_reciprocate
from compas_3gs.algorithms import volmesh_ud
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.utilities import get_force_colors_uv

OFFSET = 2
KMAX = 500
SCALE = 0.5

# ------------------------------------------------------------------------------
#   1. force diagram from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_honeycomb.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

# ------------------------------------------------------------------------------
#   2. form diagram: dual network of the force diagram, moved sideways
# ------------------------------------------------------------------------------
formdiagram = volmesh_dual_network(forcediagram, cls=FormNetwork)
formdiagram.dual = forcediagram
forcediagram.primal = formdiagram

box = formdiagram.aabb()
width = box.xsize
for node in formdiagram.nodes():
    x = formdiagram.node_attribute(node, "x")
    formdiagram.node_attribute(node, "x", x + width * OFFSET)

# ------------------------------------------------------------------------------
#   3. reciprocate: form edges perpendicular to the force faces
# ------------------------------------------------------------------------------
formdiagram.update_angle_deviations()
initial_deviation = max(formdiagram.edges_attribute("_a"))

volmesh_reciprocate(forcediagram, formdiagram, kmax=KMAX, weight=1, edge_min=0.5, edge_max=20, tolerance=0.01)

formdiagram.update_angle_deviations()
print("max angle deviation:", round(initial_deviation, 2), "->", round(max(formdiagram.edges_attribute("_a")), 3), "deg")

# ------------------------------------------------------------------------------
#   4. unified diagram
# ------------------------------------------------------------------------------
uv_c_dict = get_force_colors_uv(forcediagram, formdiagram, gradient=True)

cells, prisms = volmesh_ud(forcediagram, formdiagram, scale=SCALE)

# ------------------------------------------------------------------------------
#   5. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# grey background, so that the white edge of the smallest force stays visible
viewer.config.renderer.backgroundcolor = Color.grey().lightened(50)

for cell in cells:
    mesh = Mesh.from_vertices_and_faces(cells[cell]["vertices"], cells[cell]["faces"])
    viewer.scene.add(mesh, facecolor=Color.black(), show_lines=False, name="cell {}".format(cell))

for edge in prisms:
    mesh = Mesh.from_vertices_and_faces(prisms[edge]["vertices"], prisms[edge]["faces"])
    viewer.scene.add(mesh, facecolor=Color.from_rgb255(*uv_c_dict[edge]), show_lines=False, name="prism {}".format(edge))

# form diagram: the edges in the colours of their forces, as lines (a graph object draws all its edges in one colour)
for (u, v), color in uv_c_dict.items():
    line = Line(formdiagram.node_coordinates(u), formdiagram.node_coordinates(v))
    viewer.scene.add(line, linecolor=Color.from_rgb255(*color), linewidth=3, name="form edge {}-{}".format(u, v))

# camera on the unified diagram and the form diagram (the viewer opens looking at the origin)
corners = bounding_box([xyz for cell in cells for xyz in cells[cell]["vertices"]] + formdiagram.nodes_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10

viewer.show()

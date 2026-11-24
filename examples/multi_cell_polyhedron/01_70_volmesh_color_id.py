"""Colour the corresponding elements of a polyhedral form and force diagram.

The force diagram is reciprocated to its dual form diagram. Then every form edge gets an index
colour, and the interior face of the force diagram that is perpendicular to it gets the same colour.
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
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.utilities import get_force_colors_hf
from compas_3gs.utilities import get_index_colordict

OFFSET = 2
KMAX = 1000

# ------------------------------------------------------------------------------
#   1. force diagram from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_honeycomb.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

# ------------------------------------------------------------------------------
#   2. dual network (form diagram), moved sideways
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
#   3. reciprocate
# ------------------------------------------------------------------------------
formdiagram.update_angle_deviations()
initial_deviation = max(formdiagram.edges_attribute("_a"))

volmesh_reciprocate(forcediagram, formdiagram, kmax=KMAX, weight=1, edge_min=0.5, edge_max=20, tolerance=0.01, print_result_info=True)

formdiagram.update_angle_deviations()
print("max angle deviation:", round(initial_deviation, 3), "->", round(max(formdiagram.edges_attribute("_a")), 3), "degrees")

# ------------------------------------------------------------------------------
#   4. colour id: one index colour per form edge and its interior force face
# ------------------------------------------------------------------------------
uv_c_dict = get_index_colordict(list(formdiagram.edges()))
hf_c_dict = get_force_colors_hf(forcediagram, formdiagram, uv_c_dict=uv_c_dict)

faces_to_draw = [hfkey for hfkey in forcediagram.faces() if not forcediagram.is_halfface_on_boundary(hfkey)]

# ------------------------------------------------------------------------------
#   5. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# force diagram: the edges of every cell
for cell in forcediagram.cells():
    viewer.scene.add(forcediagram.cell_to_mesh(cell), show_faces=False, show_lines=True, name="cell {}".format(cell))

# force diagram: the interior faces in their colours, the boundary faces are not drawn
faces = {hfkey: forcediagram.halfface_vertices(hfkey) for hfkey in faces_to_draw}
vertices = {vkey: forcediagram.vertex_coordinates(vkey) for face in faces.values() for vkey in face}
facecolor = {hfkey: Color.from_rgb255(*hf_c_dict[hfkey]) for hfkey in faces}
viewer.scene.add(Mesh.from_vertices_and_faces(vertices, faces), facecolor=facecolor, show_lines=False, name="interior faces")

# form diagram: the edges in their colours, as lines (a graph object draws all its edges in one colour)
for (u, v), color in uv_c_dict.items():
    viewer.scene.add(Line(formdiagram.node_coordinates(u), formdiagram.node_coordinates(v)), linecolor=Color.from_rgb255(*color), linewidth=3, name="edge {}-{}".format(u, v))

# camera on both diagrams (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz") + formdiagram.nodes_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10

viewer.show()

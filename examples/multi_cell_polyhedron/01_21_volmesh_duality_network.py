"""Make the dual network of a force volmesh and show how their elements correspond.

Every cell of the volmesh is a node of the dual network, and every interior face between two cells
is an edge. The directed edge (u, v) and the halfface of cell u that points to cell v are drawn in
the same colour, with the halfface normal as an arrow along the edge direction.
"""

import json
import os

from compas.colors import Color
from compas.geometry import Point
from compas.geometry import Polygon
from compas.geometry import Vector
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.utilities import get_index_colordict
from compas_3gs.utilities import pair_uv_to_hf

OFFSET = 3
FACE_NORMAL_SCALE = 1.0

# ------------------------------------------------------------------------------
#   1. force volmesh from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_cubes_4x4.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

# ------------------------------------------------------------------------------
#   2. dual network, moved sideways by OFFSET times its width
# ------------------------------------------------------------------------------
dual_network = volmesh_dual_network(forcediagram, cls=FormNetwork)

width = dual_network.aabb().xsize
for node in dual_network.nodes():
    x = dual_network.node_attribute(node, "x")
    dual_network.node_attribute(node, "x", x + width * OFFSET)

interior_faces = [hfkey for hfkey in forcediagram.faces() if not forcediagram.is_halfface_on_boundary(hfkey)]

print("cells:", forcediagram.number_of_cells(), "-> nodes:", dual_network.number_of_nodes())
print("interior faces:", len(interior_faces), "-> edges:", dual_network.number_of_edges())

# ------------------------------------------------------------------------------
#   3. directed edges (u, v) and the halffaces of cell u that point to cell v
# ------------------------------------------------------------------------------
uv_hf_dict = pair_uv_to_hf(dual_network, forcediagram)

for uv, hfkey in uv_hf_dict.items():
    print("edge", uv, "-> halfface", hfkey)

cell_c_dict = get_index_colordict(list(forcediagram.cells()))
uv_c_dict = get_index_colordict(list(dual_network.edges()))

# ------------------------------------------------------------------------------
#   4. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# force diagram: edges of every cell, with the cell keys as tags
for cell in forcediagram.cells():
    mesh = forcediagram.cell_to_mesh(cell)
    viewer.scene.add(mesh, show_faces=False, show_lines=True, name="cell {}".format(cell))
    viewer.scene.add(Tag(str(cell), forcediagram.cell_centroid(cell), color=Color.from_rgb255(*cell_c_dict[cell]), height=20))

# dual network: node keys as tags, in the colours of their cells
for node in dual_network.nodes():
    viewer.scene.add(Tag(str(node), dual_network.node_coordinates(node), color=Color.from_rgb255(*cell_c_dict[node]), height=20))

# directed edges as arrows, and their halffaces with the halfface normals, in the same colours
for (u, v), hfkey in uv_hf_dict.items():
    color = Color.from_rgb255(*uv_c_dict[(u, v)])
    start = dual_network.node_coordinates(u)
    edge = subtract_vectors(dual_network.node_coordinates(v), start)
    viewer.scene.add(Vector(*edge), anchor=Point(*start), linecolor=color, linewidth=3, name="edge {}-{}".format(u, v))

    center = forcediagram.face_centroid(hfkey)
    normal = scale_vector(forcediagram.halfface_normal(hfkey), FACE_NORMAL_SCALE)
    viewer.scene.add(Polygon(forcediagram.face_coordinates(hfkey)), surfacecolor=color, show_lines=False, opacity=0.5, name="halfface {}".format(hfkey))
    viewer.scene.add(Vector(*normal), anchor=Point(*center), linecolor=color, linewidth=3, name="normal {}".format(hfkey))

# camera on both diagrams (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz") + dual_network.nodes_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10

viewer.show()

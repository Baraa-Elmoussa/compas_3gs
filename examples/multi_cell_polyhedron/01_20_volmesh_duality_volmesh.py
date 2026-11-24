"""Construct the dual volmesh of a force volmesh.

Every cell of the force volmesh becomes a vertex of the dual volmesh (same key, same colour), and
every interior vertex a cell. Each dual edge (u, v) is drawn as an arrow from u to v, and the halfface
of cell u that points to cell v as a face with its normal, in the same colour as the edge.
"""

import json
import os

from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas.colors import Color
from compas.geometry import Point
from compas.geometry import Polygon
from compas.geometry import Vector
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_dual_volmesh
from compas_3gs.diagrams import ForceVolMesh
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
#   2. dual volmesh, moved sideways by OFFSET times its width
# ------------------------------------------------------------------------------
dual_volmesh = volmesh_dual_volmesh(forcediagram)

box = dual_volmesh.aabb()
width = box.xsize
for vkey in dual_volmesh.vertices():
    x = dual_volmesh.vertex_attribute(vkey, "x")
    dual_volmesh.vertex_attribute(vkey, "x", x + width * OFFSET)

print("dual volmesh: cells", dual_volmesh.number_of_cells(), "vertices", dual_volmesh.number_of_vertices(), "edges", dual_volmesh.number_of_edges())

# ------------------------------------------------------------------------------
#   3. directed halffaces of the force volmesh and directed edges of the dual volmesh
# ------------------------------------------------------------------------------
cell_c_dict = get_index_colordict(list(forcediagram.cells()))
uv_c_dict = get_index_colordict(list(dual_volmesh.edges()))

uv_hf_dict = pair_uv_to_hf(dual_volmesh, forcediagram)

# ------------------------------------------------------------------------------
#   4. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# force volmesh: edges of each cell, cell keys at the cell centroids
for cell in forcediagram.cells():
    mesh = forcediagram.cell_to_mesh(cell)
    viewer.scene.add(mesh, show_faces=False, show_lines=True, name="cell {}".format(cell))
    viewer.scene.add(Tag(str(cell), forcediagram.cell_centroid(cell), color=Color.from_rgb255(*cell_c_dict[cell]), height=20))

# dual volmesh: faces of each cell, vertex keys (= force cell keys) in the colour of their force cell
for cell in dual_volmesh.cells():
    mesh = dual_volmesh.cell_to_mesh(cell)
    viewer.scene.add(mesh, facecolor={fkey: Color.white() for fkey in mesh.faces()}, show_lines=False, opacity=0.5, name="dual cell {}".format(cell))

for vkey in dual_volmesh.vertices():
    viewer.scene.add(Tag(str(vkey), dual_volmesh.vertex_coordinates(vkey), color=Color.from_rgb255(*cell_c_dict[vkey]), height=20))

# directed dual edges as arrows from u to v, and their halffaces with the halfface normals, in the same colours
for (u, v), hfkey in uv_hf_dict.items():
    color = Color.from_rgb255(*uv_c_dict[(u, v)])
    start = dual_volmesh.vertex_coordinates(u)
    edge = subtract_vectors(dual_volmesh.vertex_coordinates(v), start)
    viewer.scene.add(Vector(*edge), anchor=Point(*start), linecolor=color, linewidth=3, name="edge {}-{}".format(u, v))

    center = forcediagram.face_centroid(hfkey)
    normal = scale_vector(forcediagram.halfface_normal(hfkey), FACE_NORMAL_SCALE)
    viewer.scene.add(Polygon(forcediagram.face_coordinates(hfkey)), surfacecolor=color, show_lines=False, opacity=0.5, name="halfface {}".format(hfkey))
    viewer.scene.add(Vector(*normal), anchor=Point(*center), linecolor=color, linewidth=3, name="normal {}".format(hfkey))

# camera on both diagrams (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz") + dual_volmesh.vertices_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10

viewer.show()

"""Reciprocate a polyhedral force diagram and its form diagram.

The form diagram is the dual network of the force diagram: one node per cell and one edge per
interior face. Reciprocation makes every form edge perpendicular to its face in the force diagram,
by moving the form diagram (WEIGHT = 1), the force diagram (WEIGHT = 0), or both (in between).
"""

import json
import os

from compas.colors import Color
from compas.geometry import angle_vectors
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas_viewer import Viewer

from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.algorithms import volmesh_reciprocate
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork

OFFSET = 2
WEIGHT = 1.0
KMAX = 500
EDGE_MIN = 0.5
EDGE_MAX = 20
TOLERANCE = 0.01

# ------------------------------------------------------------------------------
#   1. force diagram from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_honeycomb.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

# ------------------------------------------------------------------------------
#   2. form diagram: the dual network, moved sideways
# ------------------------------------------------------------------------------
formdiagram = volmesh_dual_network(forcediagram, cls=FormNetwork)

box = formdiagram.aabb()
width = box.xsize
for node in formdiagram.nodes():
    x = formdiagram.node_attribute(node, "x")
    formdiagram.node_attribute(node, "x", x + width * OFFSET)

# ------------------------------------------------------------------------------
#   3. reciprocate, recording every 2nd iteration
# ------------------------------------------------------------------------------
initial_angle = max(angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True)
                    for u, v in formdiagram.edges())

frames = [({vkey: forcediagram.vertex_coordinates(vkey) for vkey in forcediagram.vertices()},
           {node: formdiagram.node_coordinates(node) for node in formdiagram.nodes()})]


def callback(forcediagram, formdiagram, k, args):
    if k % 2 == 0:
        frames.append(({vkey: forcediagram.vertex_coordinates(vkey) for vkey in forcediagram.vertices()},
                       {node: formdiagram.node_coordinates(node) for node in formdiagram.nodes()}))


volmesh_reciprocate(forcediagram,
                    formdiagram,
                    kmax=KMAX,
                    weight=WEIGHT,
                    edge_min=EDGE_MIN,
                    edge_max=EDGE_MAX,
                    tolerance=TOLERANCE,
                    callback=callback,
                    print_result_info=True)

frames.append(({vkey: forcediagram.vertex_coordinates(vkey) for vkey in forcediagram.vertices()},
               {node: formdiagram.node_coordinates(node) for node in formdiagram.nodes()}))

final_angle = max(angle_vectors(forcediagram.halfface_normal(forcediagram.cell_pair_halffaces(u, v)[0]), formdiagram.edge_direction((u, v)), deg=True)
                  for u, v in formdiagram.edges())

print("max angle between form edges and face normals:", round(initial_angle, 4), "->", round(final_angle, 4), "degrees")

# ------------------------------------------------------------------------------
#   4. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# force diagram: one mesh per cell
cellmeshes = {}
cellobjs = {}
for cell in forcediagram.cells():
    cellmeshes[cell] = forcediagram.cell_to_mesh(cell)
    cellobjs[cell] = viewer.scene.add(cellmeshes[cell], facecolor=Color.grey().lightened(50), show_lines=True, opacity=0.5, name="force cell {}".format(cell))

# form diagram
formobj = viewer.scene.add(formdiagram, edgecolor=Color.black(), nodecolor=Color.black(), show_points=True, linewidth=2, pointsize=5, name="form diagram")

# camera on both diagrams (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz") + formdiagram.nodes_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10


@viewer.on(interval=20, frames=len(frames))
def reciprocate(frame):
    force_xyz, form_xyz = frames[frame]
    for vkey, xyz in force_xyz.items():
        forcediagram.vertex_attributes(vkey, "xyz", xyz)
    for cell, mesh in cellmeshes.items():
        for index, vkey in enumerate(forcediagram.cell_vertices(cell)):
            mesh.vertex_attributes(index, "xyz", force_xyz[vkey])
        cellobjs[cell].update(update_data=True)
    for node, xyz in form_xyz.items():
        formdiagram.node_attributes(node, "xyz", xyz)
    formobj.update(update_data=True)


viewer.show()

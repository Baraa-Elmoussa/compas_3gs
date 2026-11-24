"""Planarise the faces of a volmesh, with some of its vertices fixed.

The faces are coloured by their flatness relative to the start: red is as warped as at the
start, white is flat. Faces that start flat (all interior faces among them) stay white, even
while the run warps them a little. The fixed vertices are the black points.
"""

import json
import os

from compas.colors import Color
from compas.geometry import Point
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas_viewer import Viewer

from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_planarise
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.utilities import compare_initial_current
from compas_3gs.utilities import volmesh_face_flatness

KMAX = 500
TOLERANCE_FLAT = 0.05

# the six lowest vertices of the diagram (z = -7.125), as if it stood on them
FIX_VKEYS = [0, 1, 12, 13, 20, 30]

# ------------------------------------------------------------------------------
#   1. force diagram from data: 9 cells with warped faces
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_warped.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

# ------------------------------------------------------------------------------
#   2. planarise, recording every 10th iteration
# ------------------------------------------------------------------------------
initial_flatness = volmesh_face_flatness(forcediagram)
fixed_xyz = {vkey: forcediagram.vertex_coordinates(vkey) for vkey in FIX_VKEYS}

frames = []


def callback(forcediagram, k, args):
    if k % 10:
        return
    current_flatness = volmesh_face_flatness(forcediagram)
    colors = compare_initial_current(current_flatness,
                                     initial_flatness,
                                     color_scheme=lambda i: Color(1.0, 1 - min(1.0, i), 1 - min(1.0, i)).rgb255)
    # the colours are keyed by one halfface per face: give the opposite halfface the same colour
    facecolor = {}
    for hfkey, color in colors.items():
        facecolor[hfkey] = color
        if not forcediagram.is_halfface_on_boundary(hfkey):
            facecolor[forcediagram.halfface_opposite_halfface(hfkey)] = color
    xyz = {vkey: forcediagram.vertex_coordinates(vkey) for vkey in forcediagram.vertices()}
    frames.append((xyz, facecolor))


# the initial state
callback(forcediagram, 0, None)

volmesh_planarise(forcediagram,
                  kmax=KMAX,
                  fix_vkeys=FIX_VKEYS,
                  fix_boundary_normals=False,
                  tolerance_flat=TOLERANCE_FLAT,
                  callback=callback,
                  print_result_info=True)

# the final state (the algorithm does not call back at the iteration that converges)
callback(forcediagram, 0, None)

# the tolerance is checked against the plane through the face centre along the face normal,
# the flatness printed here against the best-fit plane, so the two can differ a little
print("max flatness:", round(max(initial_flatness.values()), 4), "->", round(max(volmesh_face_flatness(forcediagram).values()), 4))
print("largest move of a fixed vertex:", max(distance_point_point(fixed_xyz[vkey], forcediagram.vertex_coordinates(vkey)) for vkey in FIX_VKEYS))

# ------------------------------------------------------------------------------
#   3. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# one mesh per cell: mesh face i is cell_faces(cell)[i], mesh vertex j is cell_vertices(cell)[j]
cellobjs = {}
for cell in forcediagram.cells():
    mesh = forcediagram.cell_to_mesh(cell)
    cellobjs[cell] = viewer.scene.add(mesh, facecolor={fkey: Color.white() for fkey in mesh.faces()}, show_lines=True, opacity=0.8, name="cell {}".format(cell))

for vkey in FIX_VKEYS:
    viewer.scene.add(Point(*forcediagram.vertex_coordinates(vkey)), pointcolor=Color.black(), pointsize=15)

# camera on the diagram (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10


@viewer.on(interval=250, frames=len(frames))
def planarise(frame):
    xyz, facecolor = frames[frame]
    for cell, obj in cellobjs.items():
        for j, vkey in enumerate(forcediagram.cell_vertices(cell)):
            obj.mesh.vertex_attributes(j, "xyz", xyz[vkey])
        obj.facecolor = {i: facecolor.get(hfkey, Color.white()) for i, hfkey in enumerate(forcediagram.cell_faces(cell))}
        obj.update(update_data=True)


viewer.show()

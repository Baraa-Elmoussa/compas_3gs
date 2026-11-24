"""Arearise chosen halffaces of a volmesh to a target area.

The top and bottom halffaces of the centre cell of a 3x3 mat are resized to a target area while
all face normals are kept. The outer boundary does not move, so the other cells change size to make
room (the corner cells shrink, the edge cells get longer). The target halffaces are coloured by
their area deviation relative to the start: blue is as far from the target as at the start, white
is at the target.
"""

import json
import os

from compas.colors import Color
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_planarise
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.utilities import compare_initial_current
from compas_3gs.utilities import volmesh_face_areaness

HFKEYS = [16, 17]  # the bottom and top halffaces of the centre cell (cell 2)
TARGET_AREA = 200.0
KMAX = 1000

# ------------------------------------------------------------------------------
#   1. force diagram from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_mat.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

# ------------------------------------------------------------------------------
#   2. target areas of the chosen halffaces
# ------------------------------------------------------------------------------
target_areas = {hfkey: TARGET_AREA for hfkey in HFKEYS}

initial_areas = {hfkey: forcediagram.halfface_area(hfkey) for hfkey in HFKEYS}
initial_areaness = volmesh_face_areaness(forcediagram, target_areas)

# ------------------------------------------------------------------------------
#   3. arearise with all normals fixed, recording every 10th iteration
# ------------------------------------------------------------------------------
frames = []


def callback(forcediagram, k, args):
    if k % 10 == 0:
        current_areaness = volmesh_face_areaness(forcediagram, target_areas)
        facecolor = compare_initial_current(current_areaness,
                                            initial_areaness,
                                            color_scheme=lambda i: Color(1 - min(1.0, i), 1 - min(1.0, i), 1.0).rgb255)
        xyz = {vkey: forcediagram.vertex_coordinates(vkey) for vkey in forcediagram.vertices()}
        frames.append((xyz, facecolor))


# the initial state
callback(forcediagram, 0, None)

volmesh_planarise(forcediagram,
                  kmax=KMAX,
                  target_areas=target_areas,
                  fix_all_normals=True,
                  tolerance_area=0.01,
                  callback=callback,
                  print_result_info=True)

# the final state (the algorithm does not call back at the iteration that converges)
callback(forcediagram, 0, None)

final_facecolor = frames[-1][1]

for hfkey in HFKEYS:
    print("halfface", hfkey, "target:", TARGET_AREA, "area:", round(initial_areas[hfkey], 3), "->", round(forcediagram.halfface_area(hfkey), 3))

# ------------------------------------------------------------------------------
#   4. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# one mesh per cell: mesh face i is halfface cell_faces(cell)[i], mesh vertex j is vertex cell_vertices(cell)[j]
meshes = {}
cellobjs = {}
for cell in forcediagram.cells():
    meshes[cell] = forcediagram.cell_to_mesh(cell)
    facecolor = {index: final_facecolor.get(hfkey, Color.white()) for index, hfkey in enumerate(forcediagram.cell_faces(cell))}
    cellobjs[cell] = viewer.scene.add(meshes[cell], facecolor=facecolor, show_lines=True, name="cell {}".format(cell))

# initial and final areas of the target halffaces
for hfkey in HFKEYS:
    text = "{} -> {}".format(round(initial_areas[hfkey], 2), round(forcediagram.halfface_area(hfkey), 2))
    viewer.scene.add(Tag(text, forcediagram.face_centroid(hfkey), color=Color.black(), height=20))

# camera on the diagram (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10


@viewer.on(interval=250, frames=len(frames))
def arearise(frame):
    xyz, facecolor = frames[frame]
    for vkey in xyz:
        forcediagram.vertex_attributes(vkey, "xyz", xyz[vkey])
    for cell, cellobj in cellobjs.items():
        for index, vkey in enumerate(forcediagram.cell_vertices(cell)):
            meshes[cell].vertex_attributes(index, "xyz", forcediagram.vertex_coordinates(vkey))
        cellobj.facecolor = {index: facecolor.get(hfkey, Color.white()) for index, hfkey in enumerate(forcediagram.cell_faces(cell))}
        cellobj.update(update_data=True)


viewer.show()

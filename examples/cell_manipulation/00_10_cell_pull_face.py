"""Pull a face of a cell along its normal.

The face moves to planes parallel to itself; its vertices slide along the neighbouring edges.
"""

import json
import os

from compas.colors import Color
from compas.geometry import add_vectors
from compas.geometry import scale_vector
from compas_viewer import Viewer

from compas_3gs import DATA
from compas_3gs.diagrams import Cell
from compas_3gs.operations import cell_relocate_face

FKEY = 5
DISTANCE = -1.0
STEPS = 50

# ------------------------------------------------------------------------------
#   1. cell from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "cells.json")) as f:
    data = json.load(f)["cells"][1]

cell = Cell.from_vertices_and_faces(data["vertices"], data["faces"])

# ------------------------------------------------------------------------------
#   2. pull the face, one step per frame
# ------------------------------------------------------------------------------
center = cell.face_centroid(FKEY)
normal = cell.face_normal(FKEY)
area = cell.face_area(FKEY)

frames = []
for i in range(STEPS + 1):
    xyz = add_vectors(center, scale_vector(normal, DISTANCE * i / STEPS))
    cell_relocate_face(cell, FKEY, xyz, normal)
    frames.append({vkey: cell.vertex_coordinates(vkey) for vkey in cell.vertices()})

print("face", FKEY, "area:", round(area, 3), "->", round(cell.face_area(FKEY), 3))

# ------------------------------------------------------------------------------
#   3. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

facecolor = {fkey: Color.white() for fkey in cell.faces()}
facecolor[FKEY] = Color.from_hex("#0092d2")
cellobj = viewer.scene.add(cell, facecolor=facecolor, show_lines=True, name="cell")


@viewer.on(interval=40, frames=len(frames))
def pull(frame):
    for vkey, xyz in frames[frame].items():
        cell.vertex_attributes(vkey, "xyz", xyz)
    cellobj.update(update_data=True)


viewer.show()

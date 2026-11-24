"""Arearise one face of a cell: move it along its normal until it reaches a target area."""

import json
import os

from compas_viewer import Viewer

from compas.colors import Color
from compas.geometry import dot_vectors
from compas_3gs import DATA
from compas_3gs.algorithms import cell_arearise_face
from compas_3gs.diagrams import Cell
from compas_3gs.operations import cell_relocate_face

FKEY = 0
TARGET_AREA = 5.0

# ------------------------------------------------------------------------------
#   1. cell from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "cell_gss.json")) as f:
    data = json.load(f)

cell = Cell.from_vertices_and_faces(data["vertices"], data["faces"])

area = cell.face_area(FKEY)
center = cell.face_centroid(FKEY)
normal = cell.face_normal(FKEY)

# ------------------------------------------------------------------------------
#   2. arearise the face, recording a frame at every evaluation
# ------------------------------------------------------------------------------
frames = []


def callback(cell, args):
    i = min(1.0, abs(cell.face_area(FKEY) - TARGET_AREA) / TARGET_AREA)
    xyz = {vkey: cell.vertex_coordinates(vkey) for vkey in cell.vertices()}
    frames.append((xyz, Color(1 - i, 1 - i, 1.0)))


cell_arearise_face(cell, FKEY, TARGET_AREA, callback=callback)

# ------------------------------------------------------------------------------
#   3. check the result
# ------------------------------------------------------------------------------
new_area = cell.face_area(FKEY)
if dot_vectors(normal, cell.face_normal(FKEY)) < 0:
    new_area *= -1

if abs(new_area - TARGET_AREA) > 1:
    print("Arearisation attempted, but did not converge...")
    print("It is likely that the target area is not valid / inexistent...")
    cell_relocate_face(cell, FKEY, center, normal)
else:
    print("face", FKEY, "area:", round(area, 3), "->", round(new_area, 3), "target:", TARGET_AREA)

frames.append(({vkey: cell.vertex_coordinates(vkey) for vkey in cell.vertices()}, Color(0.0, 0.0, 1.0)))

# ------------------------------------------------------------------------------
#   4. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

cellobj = viewer.scene.add(cell, facecolor={fkey: Color.white() for fkey in cell.faces()}, show_lines=True, name="cell")


@viewer.on(interval=100, frames=len(frames))
def arearise(frame):
    xyz, color = frames[frame]
    for vkey in xyz:
        cell.vertex_attributes(vkey, "xyz", xyz[vkey])
    cellobj.facecolor = {fkey: color if fkey == FKEY else Color.white() for fkey in cell.faces()}
    cellobj.update(update_data=True)


viewer.show()

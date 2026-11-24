"""Arearise several faces of a cell at once, optionally keeping all face normals fixed.

With fixed normals the targets must be compatible with a closed cell (the area-weighted normals
sum to zero); the targets below are reached with free normals.

The faces are coloured by their flatness relative to the start: blue is as warped as at the
start, white is flat.
"""

import json
import os

from compas.colors import Color
from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas_3gs import DATA
from compas_3gs.algorithms import cell_planarise
from compas_3gs.diagrams import Cell
from compas_3gs.utilities import cell_face_flatness
from compas_3gs.utilities import compare_initial_current

TARGET_AREAS = {2: 50.0, 3: 70.0, 6: 15.0}
FIX_NORMALS = False
KMAX = 500

# ------------------------------------------------------------------------------
#   1. cell from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "cells.json")) as f:
    data = json.load(f)["cells"][0]

cell = Cell.from_vertices_and_faces(data["vertices"], data["faces"])

initial_areas = {fkey: cell.face_area(fkey) for fkey in TARGET_AREAS}

target_normals = {}
if FIX_NORMALS:
    for fkey in cell.faces():
        target_normals[fkey] = cell.face_normal(fkey)

# ------------------------------------------------------------------------------
#   2. arearise, recording a frame at every iteration
# ------------------------------------------------------------------------------
initial_flatness = cell_face_flatness(cell)

frames = []


def callback(cell, k, args):
    current_flatness = cell_face_flatness(cell)
    facecolor = compare_initial_current(current_flatness,
                                        initial_flatness,
                                        color_scheme=lambda i: Color(1 - min(1.0, i), 1 - min(1.0, i), 1.0).rgb255)
    xyz = {vkey: cell.vertex_coordinates(vkey) for vkey in cell.vertices()}
    frames.append((xyz, facecolor))


cell_planarise(cell,
               kmax=KMAX,
               target_normals=target_normals,
               target_areas=TARGET_AREAS,
               callback=callback,
               print_result_info=True)

for fkey, target in TARGET_AREAS.items():
    print("face", fkey, "area:", round(initial_areas[fkey], 3), "->", round(cell.face_area(fkey), 3), "target:", target)

# ------------------------------------------------------------------------------
#   3. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

cellobj = viewer.scene.add(cell, facecolor={fkey: Color.white() for fkey in cell.faces()}, show_lines=True, name="cell")
for fkey, target in TARGET_AREAS.items():
    viewer.scene.add(Tag("target: {}".format(target), cell.face_centroid(fkey), color=Color.black(), height=20))


@viewer.on(interval=20, frames=len(frames))
def arearise(frame):
    xyz, facecolor = frames[frame]
    for vkey in xyz:
        cell.vertex_attributes(vkey, "xyz", xyz[vkey])
    cellobj.facecolor = {fkey: facecolor.get(fkey, (255, 255, 255)) for fkey in cell.faces()}
    cellobj.update(update_data=True)


viewer.show()

"""Planarise the faces of a cell.

The faces are coloured by their flatness relative to the start: red is as warped as at the
start, white is flat.
"""

import json
import os
import random

from compas.colors import Color
from compas.geometry import add_vectors
from compas_viewer import Viewer

from compas_3gs import DATA
from compas_3gs.algorithms import cell_planarise
from compas_3gs.diagrams import Cell
from compas_3gs.utilities import cell_face_flatness
from compas_3gs.utilities import compare_initial_current

KMAX = 1000
WARP = 0.3

# ------------------------------------------------------------------------------
#   1. cell from data, with its vertices moved so that the faces are warped
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "cells.json")) as f:
    data = json.load(f)["cells"][3]

cell = Cell.from_vertices_and_faces(data["vertices"], data["faces"])

random.seed(0)
for vkey in cell.vertices():
    move = [random.uniform(-WARP, WARP) for _ in range(3)]
    cell.vertex_attributes(vkey, "xyz", add_vectors(cell.vertex_coordinates(vkey), move))

# ------------------------------------------------------------------------------
#   2. planarise, recording a frame at every iteration
# ------------------------------------------------------------------------------
initial_flatness = cell_face_flatness(cell)

frames = []


def callback(cell, k, args):
    current_flatness = cell_face_flatness(cell)
    facecolor = compare_initial_current(current_flatness,
                                        initial_flatness,
                                        color_scheme=lambda i: Color(1.0, 1 - min(1.0, i), 1 - min(1.0, i)).rgb255)
    xyz = {vkey: cell.vertex_coordinates(vkey) for vkey in cell.vertices()}
    frames.append((xyz, facecolor))


cell_planarise(cell,
               kmax=KMAX,
               callback=callback,
               print_result_info=True)

print("max flatness:", round(max(initial_flatness.values()), 4), "->", round(max(cell_face_flatness(cell).values()), 6))

# ------------------------------------------------------------------------------
#   3. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

cellobj = viewer.scene.add(cell, facecolor={fkey: Color.white() for fkey in cell.faces()}, show_lines=True, name="cell")


@viewer.on(interval=20, frames=len(frames))
def planarise(frame):
    xyz, facecolor = frames[frame]
    for vkey in xyz:
        cell.vertex_attributes(vkey, "xyz", xyz[vkey])
    cellobj.facecolor = {fkey: facecolor.get(fkey, (255, 255, 255)) for fkey in cell.faces()}
    cellobj.update(update_data=True)


viewer.show()

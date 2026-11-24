"""Construct a cell from force vectors in equilibrium.

The vectors give an EGI, the EGI a unit cell, and the cell is planarised to face areas equal
to the vector lengths and face normals equal to the vector directions. Faces that come from
crossings of the EGI arcs (zero faces) get a target area of 0.
"""

import json
import os

from compas.colors import Color
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import add_vectors
from compas.geometry import length_vector
from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas_3gs import DATA
from compas_3gs.algorithms import cell_from_egi
from compas_3gs.algorithms import cell_planarise
from compas_3gs.algorithms import egi_from_vectors
from compas_3gs.utilities import get_index_colordict

SET = 1
KMAX = 2000
COLLAPSE_EDGE_LENGTH = 0.1

# ------------------------------------------------------------------------------
#   1. force vectors from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "equilibrated_force_vectors.json")) as f:
    data = json.load(f)[SET]

origin = data["origin"]
vectors = dict(enumerate(data["vectors"]))
target_areas = {key: length_vector(vector) for key, vector in vectors.items()}

# ------------------------------------------------------------------------------
#   2. egi
# ------------------------------------------------------------------------------
egi = egi_from_vectors(vectors, origin)

# ------------------------------------------------------------------------------
#   3. unit cell
# ------------------------------------------------------------------------------
cell = cell_from_egi(egi)

# ------------------------------------------------------------------------------
#   4. arearise the cell faces, recording every 10th iteration
# ------------------------------------------------------------------------------
for fkey in cell.faces():
    if fkey not in target_areas:
        target_areas[fkey] = 0

target_normals = {fkey: egi.vertex_attribute(fkey, "normal") for fkey in cell.faces()}

frames = [{vkey: cell.vertex_coordinates(vkey) for vkey in cell.vertices()}]


def callback(cell, k, args):
    if k % 10 == 0:
        frames.append({vkey: cell.vertex_coordinates(vkey) for vkey in cell.vertices()})


cell_planarise(cell,
               kmax=KMAX,
               target_areas=target_areas,
               target_normals=target_normals,
               collapse_edge_length=COLLAPSE_EDGE_LENGTH,
               callback=callback,
               print_result_info=True)

frames.append({vkey: cell.vertex_coordinates(vkey) for vkey in cell.vertices()})

for key in vectors:
    print("face", key, "target:", round(target_areas[key], 5), "area:", round(cell.face_area(key), 5))

# ------------------------------------------------------------------------------
#   5. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

colordict = get_index_colordict(list(vectors.keys()))

# input vectors and target areas
for key, vector in vectors.items():
    end = add_vectors(origin, vector)
    viewer.scene.add(Line(origin, end), linecolor=Color.from_rgb255(*colordict[key]), linewidth=3)
    viewer.scene.add(Tag("target: {}".format(round(target_areas[key], 3)), end, color=Color.from_rgb255(*colordict[key]), height=20))

# egi: zero vertices in red
egiobj = viewer.scene.add(egi, show_faces=False, show_lines=True, name="egi")
for vkey in egi.vertices():
    if egi.vertex_attribute(vkey, "type") == "zero":
        viewer.scene.add(Point(*egi.vertex_coordinates(vkey)), pointcolor=Color.red(), pointsize=10)

# cell: zero faces in red, final face areas as tags
facecolor = {fkey: Color.from_rgb255(*colordict[fkey]) if fkey in vectors else Color.red() for fkey in cell.faces()}
cellobj = viewer.scene.add(cell, facecolor=facecolor, show_lines=True, opacity=0.8, name="cell")
for key in vectors:
    viewer.scene.add(Tag(str(round(cell.face_area(key), 3)), cell.face_centroid(key), color=Color.black(), height=20))


@viewer.on(interval=50, frames=len(frames))
def planarise(frame):
    for vkey, xyz in frames[frame].items():
        cell.vertex_attributes(vkey, "xyz", xyz)
    cellobj.update(update_data=True)


viewer.show()

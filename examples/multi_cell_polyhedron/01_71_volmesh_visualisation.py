"""Five ways to visualise a polyhedral force diagram and its reciprocal form diagram.

DISPLAY picks one: "id" gives every form edge and its face in the force diagram the same colour,
"compression_tension" colours them by force (blue compression, red tension, darker is larger),
"pipes" draws every form edge as a pipe with a radius proportional to the area of its face,
"external_forces" draws the resultant of the boundary faces of every boundary cell as an arrow at
its node (white to green for small to large), and "internal_forces" draws the force of every edge as an
arrow at both of its nodes, pointing at the node for compression.
"""

import json
import os

from compas_viewer import Viewer

from compas.colors import Color
from compas.datastructures import Mesh
from compas.geometry import Cylinder
from compas.geometry import Line
from compas.geometry import Point
from compas.geometry import Vector
from compas.geometry import add_vectors
from compas.geometry import bounding_box
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.geometry import length_vector
from compas.geometry import scale_vector
from compas_3gs import DATA
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.algorithms import volmesh_reciprocate
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.utilities import get_force_colors_hf
from compas_3gs.utilities import get_force_colors_uv
from compas_3gs.utilities import get_force_mags
from compas_3gs.utilities import get_index_colordict
from compas_3gs.utilities import pair_uv_to_hf
from compas_3gs.utilities import valuedict_to_colordict

DISPLAY = "compression_tension"  # "id", "compression_tension", "pipes", "external_forces" or "internal_forces"
OFFSET = 2
KMAX = 500
PIPE_SCALE = 0.03  # pipe radius per unit of face area
EXTERNAL_SCALE = 0.15  # arrow length per unit of external force
INTERNAL_SCALE = 0.08  # arrow length per unit of internal force

# ------------------------------------------------------------------------------
#   1. force diagram from data
# ------------------------------------------------------------------------------
with open(os.path.join(DATA, "volmesh_honeycomb.json")) as f:
    data = json.load(f)

forcediagram = ForceVolMesh.from_vertices_and_cells(data["vertices"], data["cells"])

# ------------------------------------------------------------------------------
#   2. form diagram: dual network of the force diagram, moved sideways
# ------------------------------------------------------------------------------
formdiagram = volmesh_dual_network(forcediagram, cls=FormNetwork)
formdiagram.dual = forcediagram
forcediagram.primal = formdiagram

box = formdiagram.aabb()
width = box.xsize
for node in formdiagram.nodes():
    x = formdiagram.node_attribute(node, "x")
    formdiagram.node_attribute(node, "x", x + width * OFFSET)

# ------------------------------------------------------------------------------
#   3. reciprocate: form edges perpendicular to the force faces
# ------------------------------------------------------------------------------
formdiagram.update_angle_deviations()
initial_deviation = max(formdiagram.edges_attribute("_a"))

volmesh_reciprocate(forcediagram, formdiagram, kmax=KMAX, weight=1, edge_min=0.5, edge_max=20, tolerance=0.01)

formdiagram.update_angle_deviations()
print("max angle deviation:", round(initial_deviation, 2), "->", round(max(formdiagram.edges_attribute("_a")), 3), "deg")

# ------------------------------------------------------------------------------
#   4. viewer
# ------------------------------------------------------------------------------
viewer = Viewer()

# grey background, so that the smallest forces, drawn white, stay visible
viewer.config.renderer.backgroundcolor = Color.grey().lightened(50)

if DISPLAY == "id":
    # force diagram: the edges of every cell, and one halfface per form edge in the colour of that edge
    uv_c_dict = get_index_colordict(list(formdiagram.edges()))
    hf_c_dict = get_force_colors_hf(forcediagram, formdiagram, uv_c_dict=uv_c_dict)

    for cell in forcediagram.cells():
        viewer.scene.add(forcediagram.cell_to_mesh(cell), show_faces=False, show_lines=True, name="force cell {}".format(cell))

    faces = {hfkey: forcediagram.halfface_vertices(hfkey) for hfkey in pair_uv_to_hf(formdiagram, forcediagram).values()}
    vertices = {vkey: forcediagram.vertex_coordinates(vkey) for face in faces.values() for vkey in face}
    facecolor = {hfkey: Color.from_rgb255(*hf_c_dict[hfkey]) for hfkey in faces}
    viewer.scene.add(Mesh.from_vertices_and_faces(vertices, faces), facecolor=facecolor, show_lines=False, name="force interior faces")

    # form diagram: the edges in their index colours, as lines (a graph object draws all its edges in one colour)
    for (u, v), color in uv_c_dict.items():
        line = Line(formdiagram.node_coordinates(u), formdiagram.node_coordinates(v))
        viewer.scene.add(line, linecolor=Color.from_rgb255(*color), linewidth=3, name="form edge {}-{}".format(u, v))

elif DISPLAY == "compression_tension":
    # force diagram: the edges of every cell, and one halfface per form edge in the colour of its force
    uv_c_dict = get_force_colors_uv(forcediagram, formdiagram, gradient=True)
    hf_c_dict = get_force_colors_hf(forcediagram, formdiagram, uv_c_dict=uv_c_dict)

    for cell in forcediagram.cells():
        viewer.scene.add(forcediagram.cell_to_mesh(cell), show_faces=False, show_lines=True, name="force cell {}".format(cell))

    faces = {hfkey: forcediagram.halfface_vertices(hfkey) for hfkey in pair_uv_to_hf(formdiagram, forcediagram).values()}
    vertices = {vkey: forcediagram.vertex_coordinates(vkey) for face in faces.values() for vkey in face}
    facecolor = {hfkey: Color.from_rgb255(*hf_c_dict[hfkey]) for hfkey in faces}
    viewer.scene.add(Mesh.from_vertices_and_faces(vertices, faces), facecolor=facecolor, show_lines=False, name="force interior faces")

    # form diagram: the edges in the colours of their forces, as lines (a graph object draws all its edges in one colour)
    for (u, v), color in uv_c_dict.items():
        line = Line(formdiagram.node_coordinates(u), formdiagram.node_coordinates(v))
        viewer.scene.add(line, linecolor=Color.from_rgb255(*color), linewidth=3, name="form edge {}-{}".format(u, v))

elif DISPLAY == "pipes":
    # force diagram: one mesh per cell
    for cell in forcediagram.cells():
        viewer.scene.add(forcediagram.cell_to_mesh(cell), facecolor=Color.white(), show_lines=True, opacity=0.5, name="force cell {}".format(cell))

    # form diagram: a pipe per edge, with a radius proportional to the area of its face in the force diagram
    uv_hf_dict = pair_uv_to_hf(formdiagram, forcediagram)
    for u, v in formdiagram.edges():
        line = Line(formdiagram.node_coordinates(u), formdiagram.node_coordinates(v))
        radius = forcediagram.halfface_area(uv_hf_dict[(u, v)]) * PIPE_SCALE
        viewer.scene.add(Cylinder.from_line_and_radius(line, radius), facecolor=Color.black(), show_lines=False, name="pipe {}-{}".format(u, v))

elif DISPLAY == "external_forces":
    # force diagram: one mesh per cell
    for cell in forcediagram.cells():
        viewer.scene.add(forcediagram.cell_to_mesh(cell), facecolor=Color.white(), show_lines=True, opacity=0.5, name="force cell {}".format(cell))

    # form diagram: nodes and edges
    viewer.scene.add(formdiagram, edgecolor=Color.black(), nodecolor=Color.black(), show_points=True, linewidth=2, pointsize=5, name="form diagram")

    # external force of every boundary cell: the resultant of the area-weighted normals of its boundary faces
    load_vectors = {}
    load_mags = {}
    for cell in forcediagram.cells_on_boundaries():
        normals = [forcediagram.halfface_normal(hfkey, unitized=False) for hfkey in forcediagram.cell_faces(cell) if forcediagram.is_halfface_on_boundary(hfkey)]
        load_vectors[cell] = [sum(axis) for axis in zip(*normals)]
        load_mags[cell] = length_vector(load_vectors[cell])

    colordict = valuedict_to_colordict(load_mags, color_scheme=lambda i: Color(1 - min(1.0, i), 1.0, 1 - min(1.0, i)).rgb255)

    # an arrow at the form node of the cell, pointing at the node if the convention is -1
    for cell, vector in load_vectors.items():
        start = formdiagram.node_coordinates(cell)
        end = add_vectors(start, scale_vector(vector, EXTERNAL_SCALE))
        if forcediagram.attributes["convention"] == -1:
            start, end = end, start
        color = Color.from_rgb255(*colordict[cell])
        viewer.scene.add(Vector.from_start_end(start, end), anchor=Point(*start), linecolor=color, linewidth=3, name="external force {}".format(cell))

elif DISPLAY == "internal_forces":
    # force diagram: one mesh per cell
    for cell in forcediagram.cells():
        viewer.scene.add(forcediagram.cell_to_mesh(cell), facecolor=Color.white(), show_lines=True, opacity=0.5, name="force cell {}".format(cell))

    # form diagram: at every node, an arrow along each of its edges with a length proportional to the force
    uv_f_dict = get_force_mags(forcediagram, formdiagram)
    uv_c_dict = get_force_colors_uv(forcediagram, formdiagram, gradient=True)

    for u in formdiagram.nodes():
        for v in formdiagram.neighbors(u):
            edge = (u, v) if (u, v) in uv_f_dict else (v, u)
            force = uv_f_dict[edge]
            start = formdiagram.node_coordinates(u)
            end = add_vectors(start, scale_vector(formdiagram.edge_direction((u, v)), abs(force) * INTERNAL_SCALE))
            if force < 0:  # compression: the arrow points at the node
                start, end = end, start
            color = Color.from_rgb255(*uv_c_dict[edge])
            viewer.scene.add(Vector.from_start_end(start, end), anchor=Point(*start), linecolor=color, linewidth=3, name="internal force {}-{}".format(u, v))

else:
    raise ValueError("unknown DISPLAY: {}".format(DISPLAY))

# camera on both diagrams (the viewer opens looking at the origin)
corners = bounding_box(forcediagram.vertices_attributes("xyz") + formdiagram.nodes_attributes("xyz"))
diagonal = distance_point_point(corners[0], corners[6])
viewer.renderer.camera.target = centroid_points(corners)
viewer.renderer.camera.distance = 1.5 * diagonal
viewer.renderer.camera.scale = diagonal / 10

viewer.show()

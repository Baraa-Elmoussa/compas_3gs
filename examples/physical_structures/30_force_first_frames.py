"""Force first: design a polyhedral force diagram, get the space frame it stands for.

This is how 3GS is meant to be used. The force diagram is a volmesh (cells sharing faces):
- each cell is a node of the structure;
- each interior face (shared by two cells) is a member: its area is the member force;
- each boundary face is an external force (load or reaction) on that node.
The form diagram (the space frame) is the dual network: one node per cell, one member per
interior face. Reciprocation makes every member perpendicular to its face.

Two force diagrams, both made with a 3GS operation (`volmesh_cell_subdivide_barycentric`):
- A. an irregular tetrahedron split into 4 cells around an inner point -> a tetrahedral frame
  (4 nodes, 6 members) under 4 external forces;
- B. a cube split into 6 pyramids -> an octahedral frame (6 nodes, 12 members) under 6 loads;
- C. the same force diagram after `point_reflection`: the face normals keep their direction while
  the cells turn inside out, so the dual frame carries the same forces in tension, with the loads
  pulling outward.

Checks: every member is perpendicular to its face, and every node is in equilibrium with the
member forces from `get_force_mags` and the loads from the boundary faces.

Run with the compas-3gs env:  python examples/physical_structures/30_force_first_frames.py

"""

from compas_viewer import Viewer
from compas_viewer.scene import Tag

from compas.colors import Color
from compas.geometry import Line
from compas.geometry import Translation
from compas.geometry import add_vectors
from compas.geometry import angle_vectors
from compas.geometry import centroid_points
from compas.geometry import length_vector
from compas.geometry import normalize_vector
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas_3gs.algorithms import volmesh_dual_network
from compas_3gs.algorithms import volmesh_reciprocate
from compas_3gs.diagrams import ForceVolMesh
from compas_3gs.diagrams import FormNetwork
from compas_3gs.operations import point_reflection
from compas_3gs.operations import volmesh_cell_subdivide_barycentric
from compas_3gs.utilities import get_force_mags

viewer = Viewer()

# ------------------------------------------------------------------------------
# A. tetrahedron -> tetrahedral frame;  B. cube -> octahedral frame
# ------------------------------------------------------------------------------

tetrahedron = ForceVolMesh.from_vertices_and_cells(
    [[0.0, 0.0, 0.0], [4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [1.5, 1.0, 3.5]],
    [[[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]]],
)
cube = ForceVolMesh.from_meshgrid(2.0, 2.0, 2.0, 1, 1, 1)
volmesh_cell_subdivide_barycentric(tetrahedron, 0)
volmesh_cell_subdivide_barycentric(cube, next(iter(cube.cells())))
reflected = ForceVolMesh.from_meshgrid(2.0, 2.0, 2.0, 1, 1, 1)
volmesh_cell_subdivide_barycentric(reflected, next(iter(reflected.cells())))
point_reflection(reflected)

cases = [
    ("A. tetrahedral frame", tetrahedron, [0, 0, 0]),
    ("B. octahedral frame", cube, [0, 12, 0]),
    ("C. octahedral frame, force diagram point-reflected", reflected, [0, 24, 0]),
]

for name, force, offset in cases:
    form = volmesh_dual_network(force, cls=FormNetwork)
    form.dual, force.primal = force, form  # volmesh_dual_network does not link them
    volmesh_reciprocate(force, form, kmax=500, weight=1.0, edge_min=0.1, edge_max=100, tolerance=0.001)
    form.update_angle_deviations()

    forces = get_force_mags(force, form)  # + tension, - compression
    loads = {}
    for hf in force.halffaces_on_boundaries():
        cell = force.halfface_cell(hf)
        load = scale_vector(force.halfface_normal(hf), force.halfface_area(hf) * force.attributes["convention"])
        loads[cell] = add_vectors(loads.get(cell, [0, 0, 0]), load)

    boundary = len(force.halffaces_on_boundaries())
    interior = (len(list(force.halffaces())) - boundary) // 2
    max_angle = max(form.edge_attribute(e, "_a") for e in form.edges())
    print("\n" + name)
    print("  force diagram: {} cells, {} interior faces, {} boundary faces".format(force.number_of_cells(), interior, boundary))
    print("  form diagram:  {} nodes, {} members; max angle member/face {:.4f} deg".format(form.number_of_nodes(), form.number_of_edges(), max_angle))
    for (u, v), n in forces.items():
        print("    member {}-{}: {:+.3f} ({})".format(u, v, n, "tension" if n > 0 else "compression"))
    for node in form.nodes():
        residual = loads[node]
        for nbr in form.neighbors(node):
            n = forces.get((node, nbr), forces.get((nbr, node)))
            direction = normalize_vector(subtract_vectors(form.node_coordinates(nbr), form.node_coordinates(node)))
            residual = add_vectors(residual, scale_vector(direction, n))  # tension pulls the node towards nbr
        load = [round(c, 3) for c in loads[node]]
        print("    node {}: load {} (|P| {:.3f}), equilibrium residual {:.1e}".format(node, load, length_vector(loads[node]), length_vector(residual)))

    # viewer: force cells (shifted left), form diagram with forces and loads
    shift = Translation.from_vector(add_vectors(offset, [-8, 0, 0]))
    for c in force.cells():
        mesh = force.cell_to_mesh(c)
        mesh.transform(shift)
        viewer.scene.add(mesh, facecolor=Color.from_i(c / max(1, force.number_of_cells() - 1)), opacity=0.5, show_lines=True)
    peak = max(abs(n) for n in forces.values())
    for (u, v), n in forces.items():
        a, b = add_vectors(form.node_coordinates(u), offset), add_vectors(form.node_coordinates(v), offset)
        color = Color.blue() if n > 0 else Color.red()
        viewer.scene.add(Line(a, b), linecolor=color, linewidth=1 + 6 * abs(n) / peak)
        viewer.scene.add(Tag("{:+.2f}".format(n), centroid_points([a, b]), color=color, height=14))
    for node, load in loads.items():
        p = add_vectors(form.node_coordinates(node), offset)
        viewer.scene.add(Line(p, add_vectors(p, scale_vector(load, 0.3))), linecolor=Color.black(), linewidth=2)
    centre = centroid_points([form.node_coordinates(n) for n in form.nodes()])
    inward = sum(1 for node, load in loads.items() if angle_vectors(load, subtract_vectors(centre, form.node_coordinates(node)), deg=True) < 90)
    print("  loads pointing towards the centre of the frame: {} of {}".format(inward, len(loads)))

viewer.show()

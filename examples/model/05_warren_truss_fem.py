"""A 3D Warren truss solved twice: with 3GS (`GSModel`, statics) and with a linear finite element model.

The truss has a triangular section: two bottom chords, one top chord, Warren diagonals from each
top node to the four bottom nodes of its bay, cross ties and lateral bracing in the bottom plane.
It carries a deck load on the bottom nodes, a wind load on the top nodes and a vehicle on one side
of the deck.

A. On three bearings (pinned, roller, slider) the truss is statically determinate: as many
   independent equilibrium equations as unknowns. 3GS gets the member forces from statics alone,
   so the finite element model must give the same forces for any member stiffness. It also gives
   what statics cannot: the displacements. 3GS also gives what the finite element model does not:
   the force polyhedron of every node, here kept only where it reaches the forces (`check`).
B. With a fourth bearing the truss is statically indeterminate (degree 1). 3GS says so and gives
   no forces; the finite element model does, and they change with the member stiffness: here with
   a ten times softer bracing in the first half of the span. (Scaling a whole member group, or any
   change under loads symmetric about the middle of the span, would change nothing here: the one
   self-stress state of this truss twists it, with the bearing reactions alternating around the
   corners.)

The finite element model is the direct stiffness method for pin-jointed bars: each member adds
EA/L * d d^T to the stiffness matrix (d its unit direction), the free degrees of freedom are solved
for the loads, and the axial force of a member is EA/L times its elongation.

Run with the compas-3gs env:  python examples/model/05_warren_truss_fem.py

"""

import numpy as np

from compas.colors import Color
from compas.geometry import Line
from compas_3gs.model import GSModel
from compas_3gs.viewer import GSViewer

BAYS = 6
BAY = 2.0  # m
WIDTH = 2.0  # m
HEIGHT = 1.5  # m
DECK = [0.0, 0.0, -10.0]  # kN per bottom node
WIND = [0.0, 2.0, 0.0]  # kN per top node
VEHICLE = [0.0, 0.0, -30.0]  # kN on one side of the deck, off the middle of the span
EA = 210e6 * 0.002  # kN: steel, 20 cm2

# ------------------------------------------------------------------------------
# 1. geometry
# ------------------------------------------------------------------------------

left = [[i * BAY, -WIDTH / 2, 0.0] for i in range(BAYS + 1)]
right = [[i * BAY, WIDTH / 2, 0.0] for i in range(BAYS + 1)]
top = [[(i + 0.5) * BAY, 0.0, HEIGHT] for i in range(BAYS)]

groups = {
    "bottom chords": [(left[i], left[i + 1]) for i in range(BAYS)] + [(right[i], right[i + 1]) for i in range(BAYS)],
    "top chord": [(top[i], top[i + 1]) for i in range(BAYS - 1)],
    "diagonals": [(top[i], point) for i in range(BAYS) for point in (left[i], left[i + 1], right[i], right[i + 1])],
    "cross ties": [(left[i], right[i]) for i in range(BAYS + 1)],
    "bracing": [(left[i], right[i + 1]) if i % 2 == 0 else (right[i], left[i + 1]) for i in range(BAYS)],
}

# ------------------------------------------------------------------------------
# 2. 3GS: supports, loads, equilibrium, force polyhedra
# ------------------------------------------------------------------------------

model = GSModel.from_lines([line for lines in groups.values() for line in lines])
model.set_support(node=model.node_at(left[0]), dofs="pinned")
model.set_support(node=model.node_at(right[0]), dofs="roller")
model.set_support(node=model.node_at(left[-1]), dofs="yz")  # a slider, free along the span
for point in left[1:-1] + right[1:-1]:
    model.add_force(node=model.node_at(point), vector=DECK)
for point in top:
    model.add_force(node=model.node_at(point), vector=WIND)
model.add_force(node=model.node_at(left[2]), vector=VEHICLE)

form = model.form
print("{} nodes, {} members, {} reaction components".format(form.number_of_nodes(), form.number_of_edges(), sum(len(form.node_attribute(n, "support")) for n in form.supports())))

result = model.solve()
print("A. three bearings: determinate {}, out of balance {:.1e} kN".format(result.determinate, result.residual))

model.derive_force_diagram(kmax=500)
for node, check in model.check().items():
    if check["area_error"] > 0.01 or check["angle"] > 1.0:
        # the EGI fixes the faces and their neighbours from the force directions alone; here the force sizes need other neighbours
        model.skipped[node] = "the cell does not reach the forces (area error {:.0%}, angle {:.0f} deg)".format(check["area_error"], check["angle"])
        del model.cells[node]
        del model.cell_faces[node]
print("force polyhedra: {} nodes, {} without".format(len(model.cells), len(model.skipped)))

# ------------------------------------------------------------------------------
# 3. linear finite elements, for the same supports and loads
# ------------------------------------------------------------------------------

nodes = list(form.nodes())
index = {node: i for i, node in enumerate(nodes)}
xyz = np.array([form.node_coordinates(node) for node in nodes])
soft = [(model.node_at(a), model.node_at(b)) for a, b in groups["bracing"][: BAYS // 2]]  # the bracing of the first half of the span
stiffness = {
    "uniform": {edge: EA for edge in form.edges()},
    "soft bracing": {edge: EA / 10 if edge in soft or edge[::-1] in soft else EA for edge in form.edges()},
}

K = {}
for case, ea in stiffness.items():
    K[case] = np.zeros((3 * len(nodes), 3 * len(nodes)))
    for u, v in form.edges():
        i, j = index[u], index[v]
        d = xyz[j] - xyz[i]
        length = np.linalg.norm(d)
        k = ea[(u, v)] / length * np.outer(d, d) / length**2
        K[case][3 * i : 3 * i + 3, 3 * i : 3 * i + 3] += k
        K[case][3 * j : 3 * j + 3, 3 * j : 3 * j + 3] += k
        K[case][3 * i : 3 * i + 3, 3 * j : 3 * j + 3] -= k
        K[case][3 * j : 3 * j + 3, 3 * i : 3 * i + 3] -= k

p = np.zeros(3 * len(nodes))
for node in nodes:
    p[3 * index[node] : 3 * index[node] + 3] = list(form.node_load(node))
fixed = [3 * index[node] + "xyz".index(axis) for node in form.supports() for axis in form.node_attribute(node, "support")]
free = [dof for dof in range(3 * len(nodes)) if dof not in fixed]

displacements = np.zeros(3 * len(nodes))
displacements[free] = np.linalg.solve(K["uniform"][np.ix_(free, free)], p[free])
reactions = K["uniform"] @ displacements - p
fem = {}
for u, v in form.edges():
    i, j = index[u], index[v]
    d = xyz[j] - xyz[i]
    length = np.linalg.norm(d)
    fem[(u, v)] = EA / length * float(np.dot(d / length, displacements[3 * j : 3 * j + 3] - displacements[3 * i : 3 * i + 3]))

# ------------------------------------------------------------------------------
# 4. comparison
# ------------------------------------------------------------------------------

print("\nmember forces (kN, + tension / - compression)   largest |N|    3GS - FEM")
for name, lines in groups.items():
    edges = [(model.node_at(a), model.node_at(b)) for a, b in lines]
    edges = [edge if form.has_edge(edge) else edge[::-1] for edge in edges]
    largest = max(edges, key=lambda edge: abs(fem[edge]))
    difference = max(abs(form.edge_force(edge) - fem[edge]) for edge in edges)
    print("  {:15s}  3GS {:+9.3f}   FEM {:+9.3f}   {:.1e}".format(name, form.edge_force(largest), fem[largest], difference))
zero = [edge for edge in form.edges() if abs(fem[edge]) < 1e-9]
print("  zero-force members: {} (at the unsupported corner)".format(len(zero)))

print("\nreactions (kN)                3GS                          FEM")
for node in form.supports():
    dofs = form.node_attribute(node, "support")
    values = [round(r, 9) + 0.0 for r in list(form.node_reaction(node)) + list(reactions[3 * index[node] : 3 * index[node] + 3])]
    print("  node {:>2} ({:>3})  [{:+8.3f} {:+8.3f} {:+8.3f}]   [{:+8.3f} {:+8.3f} {:+8.3f}]".format(node, dofs, *values))

print("\nlargest displacement (FEM only): {:.1f} mm".format(1000 * np.abs(displacements).max()))
for node, reason in model.skipped.items():
    print("no force polyhedron at node {}: {}".format(node, reason))

# ------------------------------------------------------------------------------
# 5. a fourth bearing: statically indeterminate
# ------------------------------------------------------------------------------

indeterminate = GSModel.from_lines([line for lines in groups.values() for line in lines])
for node in form.nodes():
    indeterminate.form.node_attributes(node, ["support", "px", "py", "pz"], form.node_attributes(node, ["support", "px", "py", "pz"]))
indeterminate.set_support(node=indeterminate.node_at(right[-1]), dofs="roller")
try:
    indeterminate.solve()
except ValueError as error:
    print("\nB. four bearings, 3GS: {}".format(error))

fixed = [3 * index[node] + "xyz".index(axis) for node in indeterminate.form.supports() for axis in indeterminate.form.node_attribute(node, "support")]
free = [dof for dof in range(3 * len(nodes)) if dof not in fixed]
forces = {}
for case, ea in stiffness.items():
    u_case = np.zeros(3 * len(nodes))
    u_case[free] = np.linalg.solve(K[case][np.ix_(free, free)], p[free])
    forces[case] = {}
    for u, v in form.edges():
        i, j = index[u], index[v]
        d = xyz[j] - xyz[i]
        length = np.linalg.norm(d)
        forces[case][(u, v)] = ea[(u, v)] / length * float(np.dot(d / length, u_case[3 * j : 3 * j + 3] - u_case[3 * i : 3 * i + 3]))
    vertical = {node: (K[case] @ u_case - p)[3 * index[node] + 2] for node in indeterminate.form.supports()}
    print("   FEM, {:12s}: vertical reactions {} kN".format(case, ", ".join("node {} {:+.3f}".format(node, r) for node, r in vertical.items())))
change = max(form.edges(), key=lambda edge: abs(forces["uniform"][edge] - forces["soft bracing"][edge]))
before, after = forces["uniform"][change], forces["soft bracing"][change]
print("   with the bracing of the first half ten times softer, member {} goes from {:+.3f} to {:+.3f} kN".format(change, before, after))

# ------------------------------------------------------------------------------
# 6. viewer: 3GS forces and force polyhedra; beside them, the FEM deflected shape and forces
# ------------------------------------------------------------------------------

viewer = GSViewer()
viewer.add_model(model)

scale = 0.05 * BAYS * BAY / np.abs(displacements).max()  # the largest displacement drawn as 5% of the span
shift = np.array([0.0, -3 * WIDTH, 0.0])
peak = max(abs(force) for force in fem.values())
deflected = viewer.scene.add_group(name="FEM: deflected shape (x{:.0f}) and member forces".format(scale))
for u, v in form.edges():
    a = xyz[index[u]] + scale * displacements[3 * index[u] : 3 * index[u] + 3] + shift
    b = xyz[index[v]] + scale * displacements[3 * index[v] : 3 * index[v] + 3] + shift
    color = Color.grey() if abs(fem[(u, v)]) < 1e-9 else Color.blue() if fem[(u, v)] > 0 else Color.red()
    viewer.scene.add(Line(a, b), parent=deflected, linecolor=color, linewidth=1.5 + 6 * abs(fem[(u, v)]) / peak)

viewer.show()

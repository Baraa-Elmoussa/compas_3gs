"""A 3D Warren truss solved with 3GS (`GSModel`, statics) and with a linear finite element model.

The truss has a triangular section: two bottom chords, one top chord, Warren diagonals from each
top node to the four bottom nodes of its bay, cross ties and lateral bracing in the bottom plane.
It carries a deck load on the bottom nodes, a wind load on the top nodes and a vehicle on one side
of the deck.

A. On three bearings (pinned, roller, slider) the truss is statically determinate: as many
   independent equilibrium equations as unknowns. 3GS gets the member forces from statics alone,
   so the finite element model must give the same forces for any member stiffness. It also gives
   what statics cannot: the displacements. 3GS also gives what the finite element model does not:
   the force polyhedron of every node with at least four forces.
B. With a fourth bearing the truss is statically indeterminate, with one self-stress state: every
   state in equilibrium with the loads is one of them plus any multiple of it, a line of states.
   - The member stiffness picks one point on that line, the elastic state. 3GS gives it from
     statics and stiffness (the least complementary energy), and it is the finite element state.
     Here the bracing of the first half of the span is `SOFT_BRACING` times softer.
   - Any other point is a state in equilibrium too: pick it in the input with `FOURTH_BEARING`, the
     vertical reaction at the fourth bearing, or move the slider in the viewer. The script prints
     the range in which no bearing lifts off, and builds the force polyhedra of the chosen state.
     At 0 kN the fourth bearing carries nothing, and the state is that of A.
   (Scaling a whole member group, or any change under loads symmetric about the middle of the span,
   would not move the elastic state here: the self-stress state twists the truss, with the bearing
   reactions alternating around the corners.)

The finite element model is the direct stiffness method for pin-jointed bars: each member adds
EA/L * d d^T to the stiffness matrix (d its unit direction), the free degrees of freedom are solved
for the loads, and the axial force of a member is EA/L times its elongation.

Run with the compas-3gs env:  python examples/model/05_warren_truss_fem.py

"""

import numpy as np
from compas_viewer.components import Slider

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

SOFT_BRACING = 10.0  # B: the bracing of the first half of the span is this many times softer (sets the elastic state)
FOURTH_BEARING = None  # B: the vertical reaction you pick at the fourth bearing, in kN; None for the elastic state

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
# 2. A. 3GS on three bearings: supports, loads, equilibrium, force polyhedra
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

model.derive_force_diagram()  # nodes without a polyhedron are listed in model.skipped, with the reason
print("force polyhedra: {} nodes, {} without".format(len(model.cells), len(model.skipped)))

# ------------------------------------------------------------------------------
# 3. linear finite elements: A with uniform stiffness, B with the softer bracing
# ------------------------------------------------------------------------------

nodes = list(form.nodes())
index = {node: i for i, node in enumerate(nodes)}
xyz = np.array([form.node_coordinates(node) for node in nodes])
soft = [(model.node_at(a), model.node_at(b)) for a, b in groups["bracing"][: BAYS // 2]]  # the bracing of the first half of the span
stiffness = {
    "uniform": {edge: EA for edge in form.edges()},
    "soft bracing": {edge: EA / SOFT_BRACING if edge in soft or edge[::-1] in soft else EA for edge in form.edges()},
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
# 4. A. comparison
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
# 5. B. a fourth bearing: a line of states in equilibrium, the elastic one, and the one you pick
# ------------------------------------------------------------------------------

indeterminate = GSModel(form.copy())  # the same members and loads
corner = model.node_at(right[-1])
indeterminate.set_support(node=corner, dofs="roller")
bearings = list(indeterminate.form.supports())

elastic = indeterminate.solve(stiffness=stiffness["soft bracing"])
elastic_forces = indeterminate.forces()
elastic_reactions = indeterminate.reactions()
print("\nB. four bearings: determinate {}, {} self-stress state".format(elastic.determinate, len(indeterminate.self_stresses)))
state_forces, state_reactions = indeterminate.self_stresses[0]
print("   self-stress state, vertical reactions: {}".format(", ".join("node {} {:+.3f}".format(node, state_reactions[node][2]) for node in bearings)))

fixed = [3 * index[node] + "xyz".index(axis) for node in bearings for axis in indeterminate.form.node_attribute(node, "support")]
free = [dof for dof in range(3 * len(nodes)) if dof not in fixed]
u_b = np.zeros(3 * len(nodes))
u_b[free] = np.linalg.solve(K["soft bracing"][np.ix_(free, free)], p[free])
difference = 0.0
for u, v in form.edges():
    i, j = index[u], index[v]
    d = xyz[j] - xyz[i]
    length = np.linalg.norm(d)
    force = stiffness["soft bracing"][(u, v)] / length * float(np.dot(d / length, u_b[3 * j : 3 * j + 3] - u_b[3 * i : 3 * i + 3]))
    difference = max(difference, abs(force - elastic_forces[(u, v)]))
print("   elastic state, bracing of the first half {:.0f} times softer: 3GS - FEM {:.1e} kN".format(SOFT_BRACING, difference))
print("   elastic state, vertical reactions: {}".format(", ".join("node {} {:+.3f}".format(node, elastic_reactions[node][2]) for node in bearings)))

# the states along the line: reaction(t) = elastic + t * self-stress; no bearing lifts off while every vertical reaction is >= 0
low, high = -np.inf, np.inf
for node in bearings:
    if state_reactions[node][2] > 0:
        low = max(low, -elastic_reactions[node][2] / state_reactions[node][2])
    elif state_reactions[node][2] < 0:
        high = min(high, -elastic_reactions[node][2] / state_reactions[node][2])
ends = sorted(elastic_reactions[corner][2] + t * state_reactions[corner][2] for t in (low, high))
print("   no bearing lifts off for a fourth bearing reaction from {:.2f} to {:.2f} kN".format(*ends))

if FOURTH_BEARING is not None:
    indeterminate.solve(stiffness=stiffness["soft bracing"], reactions={(corner, "z"): FOURTH_BEARING})
chosen = indeterminate.forces()
change = max(form.edges(), key=lambda edge: abs(chosen[edge] - elastic_forces[edge]))
print("   chosen state, vertical reactions: {}".format(", ".join("node {} {:+.3f}".format(node, round(indeterminate.form.node_reaction(node)[2], 9) + 0.0) for node in bearings)))
print("   largest change from the elastic state: member {} {:+.3f} -> {:+.3f} kN".format(change, elastic_forces[change], chosen[change]))

indeterminate.derive_force_diagram()
print("   force polyhedra of the chosen state: {} nodes, {} without".format(len(indeterminate.cells), len(indeterminate.skipped)))

# ------------------------------------------------------------------------------
# 6. viewer: A with 3GS and FEM; B in the chosen state, with a slider for the fourth bearing
# ------------------------------------------------------------------------------

viewer = GSViewer()
viewer.add_model(model)

scale = 0.05 * BAYS * BAY / np.abs(displacements).max()  # the largest displacement drawn as 5% of the span
shift = np.array([0.0, -3 * WIDTH, 0.0])
peak = max(abs(force) for force in fem.values())
deflected = viewer.scene.add_group(name="A, FEM: deflected shape (x{:.0f}) and member forces".format(scale))
for u, v in form.edges():
    a = xyz[index[u]] + scale * displacements[3 * index[u] : 3 * index[u] + 3] + shift
    b = xyz[index[v]] + scale * displacements[3 * index[v] : 3 * index[v] + 3] + shift
    force = fem[(u, v)]
    color = Color.grey() if abs(force) < 1e-9 else Color.red() if force > 0 else Color.blue()
    line = viewer.scene.add(Line(a, b), parent=deflected, linecolor=color, linewidth=1.5 + 6 * abs(force) / peak, name="{}-{}".format(u, v))
    fx, fy, fz = force * form.edge_direction((u, v))
    # the same info as the 3GS members (GSViewer.member_info), with the FEM force: select a member to see it
    line.attributes["data"] = {
        "type": "zero" if abs(force) < 1e-9 else "tension" if force > 0 else "compression",
        "force": round(force, 3) + 0.0,
        "fx": round(fx, 3) + 0.0,
        "fy": round(fy, 3) + 0.0,
        "fz": round(fz, 3) + 0.0,
    }

shift = [0.0, -6 * WIDTH, 0.0]
drawing = viewer.add_form(indeterminate.form, labels=False, offset=shift)
drawing["form"].name = "B, four bearings: the chosen state"
drawing["reactions"].show = False  # the arrows keep the state they were drawn in, while the slider changes the reactions
viewer.add_cells(indeterminate, offset=[1.5 * BAYS * BAY + 2.0, shift[1], 0.0])


def pick(slider, value):
    indeterminate.solve(stiffness=stiffness["soft bracing"], reactions={(corner, "z"): value})
    viewer.update_form(indeterminate.form, drawing)


viewer.ui.sidedock.show = True
start = round(elastic_reactions[corner][2] if FOURTH_BEARING is None else FOURTH_BEARING, 1)
viewer.ui.sidedock.add(Slider(title="Bearing 4 (kN)", value=start, min_val=round(ends[0], 1), max_val=round(ends[1], 1), step=0.1, action=pick))

viewer.show()

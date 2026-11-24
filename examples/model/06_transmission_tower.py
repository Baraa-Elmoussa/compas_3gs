"""A lattice transmission tower with two cross-arms, under conductor and wind loads.

The tower tapers from a 6 m square base to a 2 m square top over 24 m, in six panels: four legs,
horizontal struts at every level and X-bracing on every face. Two cross-arms carry the conductors.
The members are steel angles: the legs larger than the bracing and the horizontals.

Every face is X-braced, so the tower is highly statically indeterminate: statics leaves one free
parameter per self-stress state. `solve` gives the elastic state for the member sections. The
footing reactions do not depend on that choice: by statics, the two footings of each side share
W / 4 +- M / (2 b) of the weight W and of the overturning moment M about the base, b the base width.

Select a member, in the scene tree or in the view, to see its force in the info panel.

Run with the compas-3gs env:  python examples/model/06_transmission_tower.py

"""

from compas_3gs.model import GSModel
from compas_3gs.viewer import GSViewer

PANELS = 6
HEIGHT = 24.0  # m
BASE = 3.0  # m, half the base width
TOP = 1.0  # m, half the top width
ARM = 4.0  # m, the reach of a cross-arm beyond the top
E = 210e6  # kN/m2
AREAS = {"legs": 43e-4, "bracing": 12e-4, "horizontals": 9e-4, "cross-arms": 12e-4}  # m2
CONDUCTOR = [0.0, 6.0, -15.0]  # kN per cross-arm tip: the weight of the conductors and the wind on them
WIND = [0.0, 0.5, 0.0]  # kN per node of the tower

# ------------------------------------------------------------------------------
# 1. geometry: four corners per level, legs, bracing, horizontals and cross-arms
# ------------------------------------------------------------------------------

levels = []
for k in range(PANELS + 1):
    z = k * HEIGHT / PANELS
    w = BASE - (BASE - TOP) * z / HEIGHT
    levels.append([[w, w, z], [-w, w, z], [-w, -w, z], [w, -w, z]])
tips = [[TOP + ARM, 0.0, HEIGHT], [-TOP - ARM, 0.0, HEIGHT]]

groups = {"legs": [], "bracing": [], "horizontals": [], "cross-arms": []}
for k in range(PANELS):
    for i in range(4):
        j = (i + 1) % 4
        groups["legs"].append((levels[k][i], levels[k + 1][i]))
        groups["bracing"] += [(levels[k][i], levels[k + 1][j]), (levels[k][j], levels[k + 1][i])]
        groups["horizontals"].append((levels[k + 1][i], levels[k + 1][j]))
for tip, (a, b) in zip(tips, [(0, 3), (1, 2)]):
    groups["cross-arms"] += [(tip, levels[-1][a]), (tip, levels[-1][b]), (tip, levels[-2][a]), (tip, levels[-2][b])]

model = GSModel.from_lines([line for lines in groups.values() for line in lines])
form = model.form

# ------------------------------------------------------------------------------
# 2. supports, loads and the elastic state for the member sections
# ------------------------------------------------------------------------------

for point in levels[0]:
    model.set_support(node=model.node_at(point), dofs="pinned")
for tip in tips:
    model.add_force(node=model.node_at(tip), vector=CONDUCTOR)
for level in levels[1:]:
    for point in level:
        model.add_force(node=model.node_at(point), vector=WIND)

stiffness = {}
for name, lines in groups.items():
    for a, b in lines:
        edge = (model.node_at(a), model.node_at(b))
        stiffness[edge if form.has_edge(edge) else edge[::-1]] = E * AREAS[name]

result = model.solve(stiffness=stiffness)
print("{} nodes, {} members, {} self-stress states".format(form.number_of_nodes(), form.number_of_edges(), len(model.self_stresses)))
print("out of balance: {:.1e} kN".format(result.residual))

# ------------------------------------------------------------------------------
# 3. member forces per group, and the footings against the overturning moment
# ------------------------------------------------------------------------------

print("\nmember forces (kN, + tension / - compression)")
for name, lines in groups.items():
    forces = []
    for a, b in lines:
        edge = (model.node_at(a), model.node_at(b))
        forces.append(form.edge_force(edge if form.has_edge(edge) else edge[::-1]))
    print("  {:12s} {:+9.2f} to {:+9.2f}".format(name, min(forces), max(forces)))

weight = -sum(form.node_load(node)[2] for node in form.nodes())
moment = sum(form.node_load(node)[1] * form.node_coordinates(node)[2] - form.node_load(node)[2] * form.node_coordinates(node)[1] for node in form.nodes())
print("\nfootings (kN, negative: uplift): vertical reaction from 3GS, and W / 4 +- M / (2 b) with W = {:.1f} kN, M = {:.1f} kNm".format(weight, moment))
for point in levels[0]:
    node = model.node_at(point)
    hand = weight / 4 + (moment / (2 * 2 * BASE) if point[1] > 0 else -moment / (2 * 2 * BASE))
    print("  footing at ({:+.0f}, {:+.0f}): {:+8.3f}   {:+8.3f}".format(point[0], point[1], form.node_reaction(node)[2], hand))

model.derive_force_diagram(kmax=300)
print("\nforce polyhedra: {} nodes, {} without".format(len(model.cells), len(model.skipped)))

# ------------------------------------------------------------------------------
# 4. viewer: select a member to see its force
# ------------------------------------------------------------------------------

viewer = GSViewer()
viewer.add_model(model, labels=False)
viewer.show()

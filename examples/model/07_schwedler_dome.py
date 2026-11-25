"""A Schwedler dome: rings, meridians and one diagonal per panel, on a ring of pinned supports.

Schwedler domes, first built in Berlin in the 1860s, are spherical caps of straight bars: here a
30 m span and 10 m rise, 12 meridians, four rings and a crown node. The nodes carry the self-weight
of the cladding by their share of the dome surface, and in a second load case snow on one half.

Under the self-weight the dome works as a membrane: the meridians are in compression all the way
down, while the hoop forces in the rings change sign from compression near the crown to tension near
the base. Membrane theory puts that change at 51.8 degrees from the crown for a spherical dome under
its own weight. The dome is statically indeterminate; `solve` gives the elastic state for equal
member sections. Snow on one half puts the rings partly in tension and partly in compression.

Select a member, in the scene tree or in the view, to see its force in the info panel.

Run with the compas-3gs env:  python examples/model/07_schwedler_dome.py

"""

from math import acos
from math import asin
from math import cos
from math import degrees
from math import pi
from math import sin

from compas_3gs.model import GSModel
from compas_3gs.viewer import GSViewer

SPAN = 30.0  # m
RISE = 10.0  # m
MERIDIANS = 12
RINGS = 4
WEIGHT = 1.0  # kN/m2 of dome surface
SNOW = 0.8  # kN/m2 of dome surface, on the half with x > 0

# ------------------------------------------------------------------------------
# 1. geometry: rings at equal angles from the crown
# ------------------------------------------------------------------------------

radius = ((SPAN / 2) ** 2 + RISE**2) / (2 * RISE)  # of the sphere
base = asin(SPAN / 2 / radius)  # the angle of the base ring from the crown
angles = [base * (1 - k / (RINGS + 1)) for k in range(RINGS + 1)]
rings = []
for angle in angles:
    r, z = radius * sin(angle), radius * cos(angle) - (radius - RISE)
    rings.append([[r * cos(2 * pi * i / MERIDIANS), r * sin(2 * pi * i / MERIDIANS), z] for i in range(MERIDIANS)])
crown = [0.0, 0.0, RISE]

groups = {"rings": [], "meridians": [], "diagonals": []}
for k, ring in enumerate(rings):
    groups["rings"] += [(ring[i], ring[(i + 1) % MERIDIANS]) for i in range(MERIDIANS)]
    upper = rings[k + 1] if k < RINGS else None
    for i in range(MERIDIANS):
        groups["meridians"].append((ring[i], upper[i] if upper else crown))
        if upper:
            groups["diagonals"].append((ring[i], upper[(i + 1) % MERIDIANS]))

# ------------------------------------------------------------------------------
# 2. two models: self-weight, and self-weight with snow on one half
# ------------------------------------------------------------------------------

# the share of the dome surface of a node: the spherical zone halfway to the rings next to it, split between the nodes of its ring
halfway = [(a + b) / 2 for a, b in zip(angles, angles[1:] + [0.0])]  # halfway[k]: between ring k and the ring (or the crown) above
shares = [0.0] + [2 * pi * radius**2 * (cos(halfway[k]) - cos(halfway[k - 1])) / MERIDIANS for k in range(1, RINGS + 1)]
crown_share = 2 * pi * radius**2 * (1 - cos(halfway[-1]))

models = {}
for case, snow in (("self-weight", 0.0), ("half snow", SNOW)):
    model = GSModel.from_lines([line for lines in groups.values() for line in lines])
    for point in rings[0]:
        model.set_support(node=model.node_at(point), dofs="pinned")
    for k in range(1, RINGS + 1):
        for point in rings[k]:
            model.add_force(node=model.node_at(point), vector=[0.0, 0.0, -(WEIGHT + (snow if point[0] > 1e-9 else 0.0)) * shares[k]])
    model.add_force(node=model.node_at(crown), vector=[0.0, 0.0, -WEIGHT * crown_share])
    model.solve()
    models[case] = model

form = models["self-weight"].form
print("{} nodes, {} members, {} self-stress states".format(form.number_of_nodes(), form.number_of_edges(), len(models["self-weight"].self_stresses)))

# ------------------------------------------------------------------------------
# 3. ring and meridian forces against membrane theory
# ------------------------------------------------------------------------------

print("\nself-weight (kN, + tension / - compression)")
print("  ring  angle from crown  hoop force  meridian below  membrane theory hoop")
for k in range(1, RINGS + 1):
    ring = []
    for a, b in groups["rings"][k * MERIDIANS : (k + 1) * MERIDIANS]:
        edge = (models["self-weight"].node_at(a), models["self-weight"].node_at(b))
        ring.append(form.edge_force(edge if form.has_edge(edge) else edge[::-1]))
    a, b = groups["meridians"][(k - 1) * MERIDIANS]
    edge = (models["self-weight"].node_at(a), models["self-weight"].node_at(b))
    meridian = form.edge_force(edge if form.has_edge(edge) else edge[::-1])
    theory = 1 / (1 + cos(angles[k])) - cos(angles[k])  # N_theta / (a q): positive in tension
    print("  {:>4}  {:>14.1f}  {:+10.2f}  {:+14.2f}  {:>19s}".format(k, degrees(angles[k]), sum(ring) / len(ring), meridian, "tension" if theory > 0 else "compression"))
print("  membrane theory: no hoop force at {:.1f} degrees".format(degrees(acos((5**0.5 - 1) / 2))))

print("\nhalf snow (kN): largest compression and tension")
snowy = models["half snow"]
for name, lines in groups.items():
    forces = []
    for a, b in lines:
        edge = (snowy.node_at(a), snowy.node_at(b))
        forces.append(snowy.form.edge_force(edge if snowy.form.has_edge(edge) else edge[::-1]))
    print("  {:10s} {:+9.2f}  {:+9.2f}".format(name, min(forces), max(forces)))

models["self-weight"].derive_force_diagram()
print("\nforce polyhedra under the self-weight: {} nodes, {} without".format(len(models["self-weight"].cells), len(models["self-weight"].skipped)))

# ------------------------------------------------------------------------------
# 4. viewer: the self-weight case with its force polyhedra
# ------------------------------------------------------------------------------

viewer = GSViewer()
viewer.add_model(models["self-weight"], labels=False)
viewer.show()

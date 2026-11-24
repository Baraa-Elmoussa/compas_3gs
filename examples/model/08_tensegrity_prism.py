"""A tensegrity prism: three struts held apart by nine cables, in equilibrium without any load.

Kenneth Snelson's simplest tensegrity, the module of his Needle Tower: two equilateral triangles of
cables, the top one turned against the bottom one, three vertical cables and three struts.

With no load, equilibrium can only be a self-stress state: a prestress that the struts and cables
hold in balance among themselves. Statics finds it: the prism has one self-stress state, and only
when the top triangle is turned by 150 degrees; at any other angle it has none and cannot be
prestressed. The level of the prestress is free, so it is prescribed through one strut force. Then
all struts are in compression and all cables in tension, as a tensegrity needs.

Every node has four forces, one strut and three cables, so its force polyhedron is a tetrahedron.

The prism is also a mechanism to first order: it carries an external load through the stiffness its
prestress gives it, which is beyond the linear statics of this example.

Run with the compas-3gs env:  python examples/model/08_tensegrity_prism.py

"""

from math import cos
from math import pi
from math import radians
from math import sin

from compas_3gs.model import GSModel
from compas_3gs.viewer import GSViewer

RADIUS = 1.0  # m, of the cable triangles
HEIGHT = 1.5  # m
TWIST = 150  # degrees, of the top triangle against the bottom one
PRESTRESS = -10.0  # kN, the force in each strut

# ------------------------------------------------------------------------------
# 1. the twist at which the prism can be prestressed
# ------------------------------------------------------------------------------

print("twist  self-stress states")
for twist in (120, 130, 140, 150, 160, 170):
    bottom = [[RADIUS * cos(2 * pi * i / 3), RADIUS * sin(2 * pi * i / 3), 0.0] for i in range(3)]
    top = [[RADIUS * cos(2 * pi * i / 3 + radians(twist)), RADIUS * sin(2 * pi * i / 3 + radians(twist)), HEIGHT] for i in range(3)]
    struts = [(bottom[i], top[i]) for i in range(3)]
    cables = [(bottom[i], bottom[(i + 1) % 3]) for i in range(3)] + [(top[i], top[(i + 1) % 3]) for i in range(3)] + [(bottom[i], top[i - 1]) for i in range(3)]
    trial = GSModel.from_lines(struts + cables)
    trial.set_support(node=trial.node_at(bottom[0]), dofs="pinned")
    trial.set_support(node=trial.node_at(bottom[1]), dofs="yz")
    trial.set_support(node=trial.node_at(bottom[2]), dofs="z")
    trial.solve()
    print("  {:>3}  {}".format(twist, len(trial.self_stresses)))

# ------------------------------------------------------------------------------
# 2. the prism at its twist, with a prescribed prestress
# ------------------------------------------------------------------------------

bottom = [[RADIUS * cos(2 * pi * i / 3), RADIUS * sin(2 * pi * i / 3), 0.0] for i in range(3)]
top = [[RADIUS * cos(2 * pi * i / 3 + radians(TWIST)), RADIUS * sin(2 * pi * i / 3 + radians(TWIST)), HEIGHT] for i in range(3)]
groups = {
    "struts": [(bottom[i], top[i]) for i in range(3)],
    "bottom cables": [(bottom[i], bottom[(i + 1) % 3]) for i in range(3)],
    "top cables": [(top[i], top[(i + 1) % 3]) for i in range(3)],
    "vertical cables": [(bottom[i], top[i - 1]) for i in range(3)],
}

model = GSModel.from_lines([line for lines in groups.values() for line in lines])
form = model.form
model.set_support(node=model.node_at(bottom[0]), dofs="pinned")  # six reaction components: the prism cannot move as a whole
model.set_support(node=model.node_at(bottom[1]), dofs="yz")
model.set_support(node=model.node_at(bottom[2]), dofs="z")

strut = (model.node_at(bottom[0]), model.node_at(top[0]))
strut = strut if form.has_edge(strut) else strut[::-1]
result = model.solve(forces={strut: PRESTRESS})
print("\nout of balance {:.1e} kN, reactions {:.1e} kN".format(result.residual, max(model.reactions()[node].length for node in form.supports())))

print("\nmember forces (kN, + tension / - compression)")
for name, lines in groups.items():
    forces = []
    for a, b in lines:
        edge = (model.node_at(a), model.node_at(b))
        forces.append(form.edge_force(edge if form.has_edge(edge) else edge[::-1]))
    print("  {:16s} {}".format(name, "  ".join("{:+.3f}".format(force) for force in forces)))

# ------------------------------------------------------------------------------
# 3. the force polyhedra: one tetrahedron per node
# ------------------------------------------------------------------------------

model.derive_force_diagram()
for node, check in model.check().items():
    print("node {}: {} faces, area error {:.2%}, angle {:.2f} deg".format(node, model.cells[node].number_of_faces(), check["area_error"], check["angle"]))

viewer = GSViewer()
viewer.add_model(model)
viewer.show()

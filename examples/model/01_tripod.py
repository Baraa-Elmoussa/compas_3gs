"""A tripod with `GSModel`: member lines, supports and a load, solve, force polyhedra, view.

Form diagram: three legs from an apex to pinned feet, and a load at the apex.
`solve` gives the leg forces and the reactions from equilibrium; `derive_force_diagram` builds the
force polyhedron of every node that has at least four forces in three dimensions (here the apex:
three legs and the load). The feet have only a leg and a reaction, so they are skipped.

Run with the compas-3gs env:  python examples/model/01_tripod.py

"""

from math import cos
from math import radians
from math import sin

from compas_3gs.model import GSModel
from compas_3gs.viewer import GSViewer

apex = [0.0, 0.0, 3.0]
feet = [[2.0 * cos(radians(a)), 2.0 * sin(radians(a)), 0.0] for a in (0, 120, 240)]

model = GSModel.from_lines([(apex, foot) for foot in feet])
for foot in feet:
    model.set_support(node=model.node_at(foot), dofs="pinned")
model.add_force(node=model.node_at(apex), vector=[1.5, 0.0, -10.0])

model.solve()
model.derive_force_diagram()

print(model.summary())
for node, check in model.check().items():
    print("node {}: area error {:.2%}, angle {:.2f} deg".format(node, check["area_error"], check["angle"]))
for node, reason in model.skipped.items():
    print("node {}: {}".format(node, reason))

viewer = GSViewer()
viewer.add_model(model)
viewer.show()

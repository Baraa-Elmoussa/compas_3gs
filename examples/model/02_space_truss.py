"""A tetrahelix tower on a pinned, a slider and a roller support, with `GSModel`.

The base ring (three members between the supports) makes the truss determinate with six reaction
components: three at the pinned support, two at the slider (free along y) and one at the roller
(only z). Every node above the base carries the same load.

Run with the compas-3gs env:  python examples/model/02_space_truss.py

"""

from math import cos
from math import radians
from math import sin
from math import sqrt

from compas.geometry import add_vectors
from compas.geometry import centroid_points
from compas.geometry import distance_point_point
from compas.geometry import scale_vector
from compas.geometry import subtract_vectors
from compas_3gs.model import GSModel
from compas_3gs.viewer import GSViewer

LEVELS = 6

xyz = [[2.0 * cos(radians(a)), 2.0 * sin(radians(a)), 0.0] for a in (90, 210, 330)]
side = distance_point_point(xyz[0], xyz[1])
xyz.append(add_vectors(centroid_points(xyz), [0.0, 0.0, side * sqrt(2.0 / 3.0)]))
for k in range(4, 3 + LEVELS):
    # the next node of a tetrahelix: the reflection of node k - 4 through the face of the last three nodes
    xyz.append(subtract_vectors(scale_vector(add_vectors(add_vectors(xyz[k - 1], xyz[k - 2]), xyz[k - 3]), 2.0 / 3.0), xyz[k - 4]))

lines = [(xyz[0], xyz[1]), (xyz[1], xyz[2]), (xyz[2], xyz[0])]
lines += [(xyz[k], xyz[k - i]) for k in range(3, 3 + LEVELS) for i in (1, 2, 3)]

model = GSModel.from_lines(lines)
model.set_support(node=model.node_at(xyz[0]), dofs="pinned")
model.set_support(node=model.node_at(xyz[1]), dofs="xz")
model.set_support(node=model.node_at(xyz[2]), dofs="roller")
for point in xyz[3:]:
    model.add_force(node=model.node_at(point), vector=[0.5, 0.0, -2.0])

result = model.solve()
print("determinate: {}, out of balance: {:.1e}".format(result.determinate, result.residual))
print(model.summary())

model.derive_force_diagram()
for node, check in model.check().items():
    print("node {}: area error {:.2%}, angle {:.2f} deg".format(node, check["area_error"], check["angle"]))
for node, reason in model.skipped.items():
    print("node {}: {}".format(node, reason))

viewer = GSViewer()
viewer.add_model(model)
viewer.show()

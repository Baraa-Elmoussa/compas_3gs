"""What `GSModel.solve` says when statics alone cannot give the forces.

- Four legs to pinned feet: one leg too many, the truss is statically indeterminate.
- Two legs: the apex can swing about the line through the feet, a mechanism under a sideways load.
- The same two legs carry a load in their own plane: a mechanism that carries this load only.

Run with the compas-3gs env:  python examples/model/04_determinacy.py

"""

from compas_3gs.model import GSModel

apex = [0.0, 0.0, 2.0]
corners = [[1.0, 1.0, 0.0], [-1.0, 1.0, 0.0], [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0]]

model = GSModel.from_lines([(apex, corner) for corner in corners])
for corner in corners:
    model.set_support(model.node_at(corner), "pinned")
model.add_force(model.node_at(apex), [0.0, 0.0, -10.0])
try:
    model.solve()
except ValueError as error:
    print("four legs:", error)

model = GSModel.from_lines([(apex, corners[0]), (apex, corners[2])])
for corner in (corners[0], corners[2]):
    model.set_support(model.node_at(corner), "pinned")
model.add_force(model.node_at(apex), [1.0, -1.0, -10.0])
try:
    model.solve()
except ValueError as error:
    print("two legs, sideways load:", error)

model.form.node_attributes(model.node_at(apex), ["px", "py", "pz"], [0.0, 0.0, -10.0])
result = model.solve()
print("two legs, load in their plane: forces {}, determinate {}".format([round(force, 3) for force in result.forces], result.determinate))

"""What `GSModel.solve` does when statics alone has many solutions, or none.

- Four legs to pinned feet: one leg too many, the truss is statically indeterminate with one
  self-stress state. Every state in equilibrium is one of them plus any multiple of it; `solve` gives
  the elastic one for the member stiffness, or the one fixed by a prescribed force.
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
legs = list(model.form.edges())

result = model.solve()
print("four legs: determinate {}, {} self-stress state(s)".format(result.determinate, len(result.self_stresses)))
print("  self-stress state:             {}".format([round(force, 3) for force in result.self_stresses[0][0]]))
print("  equal stiffness:               {}".format([round(force, 3) for force in result.forces]))
result = model.solve(stiffness={leg: 4.0 if leg == legs[0] else 1.0 for leg in legs})
print("  first leg four times stiffer:  {}".format([round(force, 3) for force in result.forces]))
result = model.solve(forces={legs[0]: 0.0})
print("  first leg prescribed at 0 kN:  {}".format([round(force, 3) + 0.0 for force in result.forces]))

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

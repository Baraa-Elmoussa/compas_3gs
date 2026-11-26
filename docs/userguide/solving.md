# Solving

`model.solve()` finds member forces and support reactions in equilibrium with the loads, and stores
them in the form diagram: the axial force of each member in the edge attribute `_f` (tension
positive, compression negative), the reaction of each support in the node attributes `_rx`, `_ry`, `_rz`.
{ .lead }

```python
result = model.solve()
model.forces()     # {edge: force}
model.reactions()  # {node: Vector}
print(model.summary())
```

## Statics

Every node gives three equilibrium equations; the unknowns are one force per member and one reaction
per restrained direction. `solve` sets these equations up and reports what they allow
(`EquilibriumResult`):

- **Statically determinate:** one solution. Statics alone gives the forces, whatever the member
  sections (`result.determinate` is `True`).
- **Statically indeterminate:** more unknowns than independent equations. The states in equilibrium are
  one of them plus any combination of the **self-stress states**: sets of member forces and reactions in
  equilibrium without any load. `model.self_stresses` lists them, as `({edge: force}, {node: reaction})`.
- **Mechanism:** fewer independent equations than the nodes need. A mechanism that cannot carry the load
  raises `ValueError`; one that happens to carry it (a load in the plane of two legs, for example) is
  solved, with `result.determinate` `False`.

## Picking a state of an indeterminate truss

Graphic statics does not need the structure to be determinate: any state in equilibrium has its force
diagram. `solve` takes the state with the least complementary energy, the sum of `N² L / (2 EA)` over
the members, among the states that have the forces and reactions you prescribe:

```python
model.solve()                                  # nothing prescribed, equal stiffness: the elastic state of equal members
model.solve(stiffness={edge: EA})              # the elastic state for these member sections
model.solve(forces={edge: -10.0})              # a prescribed member force
model.solve(reactions={(node, "z"): 25.0})     # a prescribed reaction component
```

- With the stiffness of the real members and rigid supports, and nothing prescribed, the state is the one a
  linear finite element model gives: [Warren truss, FEM and 3GS](../examples/model_05_warren_truss_fem.md)
  shows both agreeing to 10⁻¹² kN.
- As many independent prescribed values as there are self-stress states fix the state, whatever the
  stiffness. Values that statics already fixes, or too many of them, raise `ValueError`.
- A structure in pure prestress, with no load at all, is a self-stress state: prescribe one force to set its
  level ([Tensegrity prism](../examples/model_08_tensegrity_prism.md)).

## On plain data

`model.solve` calls `truss_equilibrium` of `compas_3gs.algorithms`, which works on plain lists:
vertex coordinates, edges as index pairs, loads and supports by vertex index.

```python
from compas_3gs.algorithms import truss_equilibrium

result = truss_equilibrium(vertices, edges, loads, supports, stiffness=None, forces=None, reactions=None)
```

Next: [Force polyhedra](force_polyhedra.md).
